process ANNOTATE {
    tag "${id}"

    input:
    tuple val(id), path(assembly)

    output:
    path "${id}",           emit: dir
    path "${id}/${id}.gff", emit: gff
    tuple val("${task.process}"), val('prokka'), eval('prokka --version 2>&1 | cut -d" " -f2'), topic: versions

    script:
    """
    gzip -cdf ${assembly} > ${id}.fasta    # -f passes uncompressed files (e.g. K-12 .fna) through
    prokka \\
        --outdir ${id} \\
        --prefix ${id} \\
        --genus Escherichia --species coli \\
        --compliant \\
        --cpus ${task.cpus} \\
        ${id}.fasta
    rm ${id}.fasta
    """

    stub:
    """
    mkdir ${id}
    touch ${id}/${id}.gff
    """
}
