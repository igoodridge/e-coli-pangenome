process SUMMARISE_ANNOTATIONS {
    input:
    path pangenome      // Panaroo output directory
    path annotations    // eggNOG-mapper .emapper.annotations
    path cog_defs       // NCBI COG definitions (cog-20.def.tab)

    output:
    path 'gene_annotations.tsv', emit: genes

    script:
    """
    summarise_annotations.py \\
        --presence ${pangenome}/gene_presence_absence.csv \\
        --eggnog ${annotations} \\
        --cog-defs ${cog_defs}
    """

    stub:
    """
    touch gene_annotations.tsv
    """
}
