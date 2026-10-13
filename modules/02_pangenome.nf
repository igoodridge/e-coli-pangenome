process PANGENOME {
    input:
    path gffs, stageAs: 'gffs/*'

    output:
    path 'panaroo',                         emit: dir
    path 'panaroo/pan_genome_reference.fa', emit: ref
    tuple val("${task.process}"), val('panaroo'), eval('panaroo --version 2>&1 | cut -d" " -f2'), topic: versions

    script:
    // GFFs go in as a list file: thousands of paths on the command line can exceed the argument limit
    """
    find gffs/ -name '*.gff' | sort > gff_list.txt
    mkdir panaroo
    panaroo \\
        -i gff_list.txt \\
        -o panaroo \\
        --clean-mode moderate \\
        --remove-invalid-genes \\
        --aligner none \\
        -t ${task.cpus}
    """

    stub:
    """
    mkdir panaroo
    touch panaroo/pan_genome_reference.fa panaroo/gene_presence_absence.csv
    """
}
