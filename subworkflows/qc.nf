include { SEQKIT_STATS; DOWNLOAD_CHECKM2_DB; CHECKM2; DOWNLOAD_SKANI_DB; SKANI_SEARCH; QC_FILTER } from '../modules/00_qc'

// QC the assemblies and extra genomes, and emit those that pass as tuple(id, assembly).
// Each report (seqkit, CheckM2, skani) is generated for every genome, or, if an existing report is given,
// only for the extra genomes (which it won't contain). Databases are only needed when a tool runs,
// and are downloaded if not set.
workflow QC {
    take:
    assemblies   // channel: path(assembly)
    extras       // channel: path(<genome_id>.fasta.gz) for extra genomes, e.g. K-12
    opts         // map: pipeline params (qc_*, checkm2_db, skani_db, skani_db_url, db_cache)

    main:
    all_genomes = assemblies.mix(extras)

    seqkit_input = opts.qc_seqkit ? extras.collect() : all_genomes.collect()
    seqkit = SEQKIT_STATS(seqkit_input).report
        .mix(opts.qc_seqkit ? channel.value(opts.qc_seqkit) : channel.empty())

    checkm2_input = opts.qc_checkm2 ? extras.collect() : all_genomes.collect()
    checkm2_db = opts.checkm2_db
        ? channel.value(opts.checkm2_db)
        : DOWNLOAD_CHECKM2_DB(checkm2_input.map { _genomes -> opts.db_cache })   // only if CheckM2 runs
    checkm2 = CHECKM2(checkm2_input, checkm2_db).report
        .mix(opts.qc_checkm2 ? channel.value(opts.qc_checkm2) : channel.empty())

    skani_input = opts.qc_skani ? extras.collect() : all_genomes.collect()
    skani_db = opts.skani_db
        ? channel.value(opts.skani_db)
        : DOWNLOAD_SKANI_DB(skani_input.map { _genomes -> opts.skani_db_url }, opts.db_cache)
    skani = SKANI_SEARCH(skani_input, skani_db).report
        .mix(opts.qc_skani ? channel.value(opts.qc_skani) : channel.empty())

    extra_ids = extras.map { f -> genomeId(f) }.collect().ifEmpty([])
    QC_FILTER(seqkit.collect(), checkm2.collect(), skani.collect(), opts.qc_config, extra_ids)

    passed_ids = QC_FILTER.out.passed
        .splitCsv(header: true, sep: '\t')
        .map { row -> row.genome_id }

    passed_genomes = all_genomes
        .map { f -> tuple(genomeId(f), f) }
        .join(passed_ids)

    emit:
    assemblies = passed_genomes                                      // tuple(id, assembly) passing QC, incl. kept extras
    passed     = QC_FILTER.out.passed                                // TSV of genomes passing QC
    summary    = QC_FILTER.out.report                                // TSV of every genome with fail reasons
    generated  = SEQKIT_STATS.out.report.mix(CHECKM2.out.report, SKANI_SEARCH.out.report)   // reports made in this run
}

def genomeId(f) {
    f.name.replaceAll(/\.(fasta|fna|fa)(\.gz)?$/, '')
}
