include { SEQKIT_STATS; DOWNLOAD_CHECKM2_DB; CHECKM2; DOWNLOAD_SKANI_DB; SKANI_SEARCH; QC_FILTER } from '../modules/00_qc'

// QC the assemblies and emit those that pass as tuple(id, assembly).
// Each report (seqkit, CheckM2, skani) is used as given if its param is set, otherwise generated.
// Databases are only needed when their report is generated, and are downloaded if not set.
workflow QC {
    take:
    assemblies   // channel: path(assembly)
    opts         // map: pipeline params (qc_*, checkm2_db, skani_db, skani_db_url, db_cache)

    main:
    all_assemblies = assemblies.collect()

    seqkit = opts.qc_seqkit
        ? channel.value(opts.qc_seqkit)
        : SEQKIT_STATS(all_assemblies).report

    checkm2 = opts.qc_checkm2
        ? channel.value(opts.qc_checkm2)
        : CHECKM2(
            all_assemblies,
            opts.checkm2_db ? channel.value(opts.checkm2_db) : DOWNLOAD_CHECKM2_DB(opts.db_cache),
        ).report

    skani = opts.qc_skani
        ? channel.value(opts.qc_skani)
        : SKANI_SEARCH(
            all_assemblies,
            opts.skani_db ? channel.value(opts.skani_db) : DOWNLOAD_SKANI_DB(opts.skani_db_url, opts.db_cache),
        ).report

    QC_FILTER(seqkit, checkm2, skani, opts.qc_config)

    passed_ids = QC_FILTER.out.passed
        .splitCsv(header: true, sep: '\t')
        .map { row -> row.genome_id }

    passed_assemblies = assemblies
        .map { f -> tuple(genomeId(f), f) }
        .join(passed_ids)

    emit:
    assemblies = passed_assemblies                    // tuple(id, assembly) for genomes passing QC
    passed     = QC_FILTER.out.passed                 // TSV of genomes passing QC
    summary    = QC_FILTER.out.report                 // TSV of every genome with fail reasons
    reports    = seqkit.mix(checkm2, skani)           // the reports QC was based on
}

def genomeId(f) {
    f.name.replaceAll(/\.(fasta|fna|fa)(\.gz)?$/, '')
}
