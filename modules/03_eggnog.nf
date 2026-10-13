process DOWNLOAD_EGGNOG_DB {
    storeDir "${cache_dir}/eggnog_${db_version}"   // downloaded once, reused by later runs

    input:
    val url
    val db_version
    val cache_dir

    output:
    path 'eggnog_db'

    script:
    """
    mkdir eggnog_db
    for f in eggnog.db eggnog.db.taxids.bin eggnog.db.fieldpresence.bin \\
             eggnog.taxa.db eggnog.taxa.db.traverse.pkl eggnog_proteins.dmnd go-basic.obo; do
        curl -fL --retry 5 -C - -o eggnog_db/\$f ${url}/\$f
    done
    check_eggnog_db.py eggnog_db ${db_version}
    """

    stub:
    """
    mkdir eggnog_db
    """
}

process EGGNOG {
    input:
    path cds
    path db
    val db_version

    output:
    path 'eggnog_pangenome.emapper.*',           emit: annotations
    path 'eggnog_pangenome.emapper.annotations', emit: table
    tuple val("${task.process}"), val('eggnog-mapper'), eval('emapper.py --version 2>&1 | grep -o "emapper-[^ ]*" | sed -E "s/emapper-//; s/-[0-9a-f]{7}\$//"'), topic: versions
    tuple val("${task.process}"), val('eggnog-db'), val(db_version), topic: versions

    script:
    // --itype CDS: emapper translates the Panaroo reference (frame 1, bacterial code) before the DIAMOND search
    """
    check_eggnog_db.py ${db} ${db_version}
    emapper.py \\
        -i ${cds} \\
        --itype CDS \\
        --trans_table 11 \\
        --output eggnog_pangenome \\
        --data_dir ${db} \\
        -m diamond \\
        --cpu ${task.cpus}
    """

    stub:
    """
    touch eggnog_pangenome.emapper.annotations
    """
}