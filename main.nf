include { QC; genomeId               } from './subworkflows/qc'
include { ANNOTATE                   } from './modules/01_annotate'
include { PANGENOME                  } from './modules/02_pangenome'
include { DOWNLOAD_EGGNOG_DB; EGGNOG } from './modules/03_eggnog'

params {
    // Directory of genome assemblies (*.fasta.gz)
    input: Path

    // QC thresholds
    qc_config: Path = "${projectDir}/conf/qc.yml"

    // Existing QC reports; any left unset are generated from the assemblies
    qc_seqkit: Path?
    qc_checkm2: Path?
    qc_skani: Path?

    // QC databases, only used when generating that report; if unset they're downloaded into db_cache
    checkm2_db: Path?
    skani_db: Path?
    skani_db_url: String = 'http://faust.compbio.cs.cmu.edu/skani-files/skani_gtdb_r226-v0.3.tar.gz'

    // Reference genomes added after QC, e.g. E. coli K-12 (path or glob; .fasta/.fna/.fa, optionally .gz)
    extra_genomes: String?

    // Use only the first N assemblies (0 = all), applied before QC; extra genomes are always included
    test_n: Integer = 0

    // Existing eggNOG DB directory; if unset, the DB is downloaded into db_cache
    eggnog_db: Path?
    eggnog_db_version: String = '7.0.0'
    eggnog_db_url: String = 'https://data.cgmlab.org/eggnog-mapper/emapper-3.0/data'

    // Where downloaded databases are kept
    db_cache: String = "${projectDir}/databases"
}

workflow {
    main:
    raw_assemblies = channel.fromPath("${params.input}/*.fasta.gz", checkIfExists: true)
        .take(params.test_n ?: -1)

    QC(raw_assemblies, params)

    extra_genomes = params.extra_genomes
        ? channel.fromPath(params.extra_genomes, checkIfExists: true).map { f -> tuple(genomeId(f), f) }
        : channel.empty()

    ANNOTATE(QC.out.assemblies.mix(extra_genomes))
    PANGENOME(ANNOTATE.out.gff.collect())

    eggnog_db = params.eggnog_db
        ? channel.value(params.eggnog_db)
        : DOWNLOAD_EGGNOG_DB(params.eggnog_db_url, params.eggnog_db_version, params.db_cache)
    EGGNOG(PANGENOME.out.ref, eggnog_db, params.eggnog_db_version)

    versions = channel.topic('versions')
        .unique()
        .map { proc, tool, version -> "${proc.tokenize(':').last()}\t${tool}\t${version}" }
        .collectFile(name: 'software_versions.tsv', seed: "process\ttool\tversion", newLine: true, sort: true)

    publish:
    qc         = QC.out.passed.mix(QC.out.summary)
    qc_reports = QC.out.reports
    annotation = ANNOTATE.out.dir
    pangenome  = PANGENOME.out.dir
    eggnog     = EGGNOG.out.annotations
    versions   = versions
}

output {
    qc         { path 'qc' }
    qc_reports { path 'qc/reports' }
    annotation { path 'prokka' }
    pangenome  { path '.' }
    eggnog     { path 'eggnog' }
    versions   { path 'pipeline_info' }
}
