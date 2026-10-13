process TRANSLATE {
    input:
    path ref

    output:
    path 'pan_genome_reference.faa', emit: faa
    tuple val("${task.process}"), val('emboss'), eval('transeq -version 2>&1 | sed "s/EMBOSS://"'), topic: versions

    script:
    """
    transeq \\
        -sequence ${ref} \\
        -outseq pan_genome_reference.faa \\
        -frame 1 \\
        -trim Y
    """

    stub:
    """
    touch pan_genome_reference.faa
    """
}
