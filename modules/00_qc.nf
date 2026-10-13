process SEQKIT_STATS {
    input:
    path assemblies, stageAs: 'assemblies/*'

    output:
    path 'seqkit_stats.tsv', emit: report
    tuple val("${task.process}"), val('seqkit'), eval('seqkit version | sed "s/seqkit v//"'), topic: versions

    script:
    """
    seqkit stats --all --tabular --threads ${task.cpus} assemblies/* > seqkit_stats.tsv
    """

    stub:
    """
    (echo file; ls assemblies/*) > seqkit_stats.tsv
    """
}

process DOWNLOAD_CHECKM2_DB {
    storeDir "${cache_dir}/checkm2"   // downloaded once, reused by later runs

    input:
    val cache_dir

    output:
    path 'CheckM2_database/uniref100.KO.1.dmnd'

    script:
    """
    checkm2 database --download --path . --no_write_json_db
    """

    stub:
    """
    mkdir CheckM2_database
    touch CheckM2_database/uniref100.KO.1.dmnd
    """
}

process CHECKM2 {
    input:
    path assemblies, stageAs: 'assemblies/*'
    path db

    output:
    path 'checkm2/quality_report.tsv', emit: report
    tuple val("${task.process}"), val('checkm2'), eval('checkm2 --version'), topic: versions

    script:
    """
    checkm2 predict \\
        --threads ${task.cpus} \\
        --input assemblies \\
        --extension fasta.gz \\
        --output-directory checkm2 \\
        --database_path ${db} \\
        --remove_intermediates
    """

    stub:
    """
    mkdir checkm2
    touch checkm2/quality_report.tsv
    """
}

process DOWNLOAD_SKANI_DB {
    storeDir "${cache_dir}/skani/${url.tokenize('/').last() - '.tar.gz'}"   // downloaded once, reused by later runs

    input:
    val url
    val cache_dir

    output:
    path 'skani_db'

    script:
    """
    curl -fL --retry 5 -C - -o skani_db.tar.gz ${url}
    mkdir skani_db
    tar -xzf skani_db.tar.gz -C skani_db --strip-components 1
    rm skani_db.tar.gz
    """

    stub:
    """
    mkdir skani_db
    """
}

process SKANI_SEARCH {
    input:
    path assemblies, stageAs: 'assemblies/*'
    path db

    output:
    path 'skani_results.tsv', emit: report
    tuple val("${task.process}"), val('skani'), eval('skani -V | cut -d" " -f2'), topic: versions

    script:
    """
    skani search \\
        -d ${db} \\
        -o skani_results.tsv \\
        -t ${task.cpus} \\
        assemblies/*
    """

    stub:
    """
    touch skani_results.tsv
    """
}

process QC_FILTER {
    input:
    path seqkit
    path checkm2
    path skani
    path config

    output:
    path 'qc_passed.tsv', emit: passed
    path 'qc_report.tsv', emit: report
    tuple val("${task.process}"), val('pandas'), eval('python -c "import pandas; print(pandas.__version__)"'), topic: versions

    script:
    """
    qc_filter.py \\
        --seqkit ${seqkit} \\
        --checkm2 ${checkm2} \\
        --skani ${skani} \\
        --config ${config}
    """

    stub:
    // Real filter if the reports are real (it only takes seconds), otherwise pass every genome in the seqkit report
    """
    if [ -s ${seqkit} ] && [ -s ${checkm2} ] && [ -s ${skani} ]; then
        qc_filter.py --seqkit ${seqkit} --checkm2 ${checkm2} --skani ${skani} --config ${config}
    else
        (echo genome_id; tail -n +2 ${seqkit} | xargs -n1 basename | sed 's/\\.fasta\\.gz\$//') | tee qc_passed.tsv > qc_report.tsv
    fi
    """
}
