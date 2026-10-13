include { QC                         } from './subworkflows/qc'
include { PREPARE_GENOME             } from './modules/00_qc'
include { ANNOTATE                   } from './modules/01_annotate'
include { PANGENOME                  } from './modules/02_pangenome'
include { DOWNLOAD_EGGNOG_DB; EGGNOG } from './modules/03_eggnog'
include { SUMMARISE_ANNOTATIONS      } from './modules/04_summarise'
include { GENOME_SELECTION           } from './subworkflows/genome_selection'

params {
    // Directory of genome assemblies (*.fasta.gz)
    input: Path

    // Extra genomes, e.g. E. coli K-12: TSV with genome_id, fasta and phylogroup columns.
    // They go through QC (see extra_genomes in conf/qc.yml) and are never limited by test_n
    extra_genomes: Path?

    // Use only the first N assemblies (0 = all), applied before QC
    test_n: Integer = 0

    // QC thresholds
    qc_config: Path = "${projectDir}/conf/qc.yml"

    // Existing QC reports; any left unset are generated from the assemblies
    qc_seqkit: Path?
    qc_checkm2: Path?
    qc_skani: Path?

    // QC databases, only used when a QC tool runs; if unset they're downloaded into db_cache
    checkm2_db: Path?
    skani_db: Path?
    skani_db_url: String = 'http://faust.compbio.cs.cmu.edu/skani-files/skani_gtdb_r226-v0.3.tar.gz'

    // Existing eggNOG DB directory; if unset, the DB is downloaded into db_cache
    eggnog_db: Path?
    eggnog_db_version: String = '7.0.0'
    eggnog_db_url: String = 'https://data.cgmlab.org/eggnog-mapper/emapper-3.0/data'

    // Existing Prokka output (<prokka_dir>/<id>/<id>.gff): genomes found there aren't re-annotated
    prokka_dir: Path?

    // NCBI COG definitions, used to describe eggNOG's COG IDs
    cog_defs: Path = "${projectDir}/assets/cog-20.def.tab"

    // Genome selection: COG gene sets for distances, picking rules, and phylogroups (TSV with guuid, phylogroup)
    cog_sets: Path = "${projectDir}/conf/cog_sets.yml"
    picking: Path = "${projectDir}/conf/picking.yml"
    phylogroups: Path?

    // Where downloaded databases are kept
    db_cache: String = "${projectDir}/databases"
}

workflow {
    main:
    raw_assemblies = channel.fromPath("${params.input}/*.fasta.gz", checkIfExists: true)
        .take(params.test_n ?: -1)

    extra_genomes = params.extra_genomes
        ? channel.fromPath(params.extra_genomes)
            .splitCsv(header: true, sep: '\t')
            .map { row -> tuple(row.genome_id, file(row.fasta, checkIfExists: true)) }
        : channel.empty()
    PREPARE_GENOME(extra_genomes)

    QC(raw_assemblies, PREPARE_GENOME.out, params)

    // Reuse existing Prokka annotations where there are any; annotate the rest
    to_annotate = QC.out.assemblies.branch { id, _assembly ->
        existing: params.prokka_dir && params.prokka_dir.resolve("${id}/${id}.gff").exists()
        missing: true
    }
    existing_gffs = to_annotate.existing.map { id, _assembly -> params.prokka_dir.resolve("${id}/${id}.gff") }
    ANNOTATE(to_annotate.missing)
    PANGENOME(ANNOTATE.out.gff.mix(existing_gffs).collect())

    eggnog_db = params.eggnog_db
        ? channel.value(params.eggnog_db)
        : DOWNLOAD_EGGNOG_DB(params.eggnog_db_url, params.eggnog_db_version, params.db_cache)
    EGGNOG(PANGENOME.out.ref, eggnog_db, params.eggnog_db_version)

    SUMMARISE_ANNOTATIONS(PANGENOME.out.dir, EGGNOG.out.table, params.cog_defs)

    GENOME_SELECTION(PANGENOME.out.dir, SUMMARISE_ANNOTATIONS.out.genes, QC.out.summary, params)

    versions = channel.topic('versions')
        .unique()
        .map { proc, tool, version -> "${proc.tokenize(':').last()}\t${tool}\t${version}" }
        .collectFile(name: 'software_versions.tsv', seed: "process\ttool\tversion", newLine: true, sort: true)

    publish:
    qc           = QC.out.passed.mix(QC.out.summary)
    qc_reports   = QC.out.generated
    annotation   = ANNOTATE.out.dir
    pangenome    = PANGENOME.out.dir
    eggnog       = EGGNOG.out.annotations
    summary      = SUMMARISE_ANNOTATIONS.out.genes
    distances    = GENOME_SELECTION.out.distances
    picks        = GENOME_SELECTION.out.picks
    plots        = GENOME_SELECTION.out.plots
    picks_summary = GENOME_SELECTION.out.summary
    versions     = versions
}

output {
    qc           { path 'qc' }
    qc_reports   { path 'qc/reports' }
    annotation   { path 'prokka' }
    pangenome    { path '.' }
    eggnog       { path 'eggnog' }
    summary      { path 'annotations' }
    distances    { path 'genome_selection' }
    picks        { path 'genome_selection' }
    plots        { path 'genome_selection' }
    picks_summary { path 'genome_selection' }
    versions     { path 'pipeline_info' }
}
