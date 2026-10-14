process GENOME_DISTANCES {
    input:
    path pangenome      // Panaroo output directory
    path annotations    // gene_annotations.tsv
    path cog_sets       // COG sets YAML

    output:
    path 'distances', emit: distances

    script:
    """
    genome_distances.py \\
        --presence ${pangenome}/gene_presence_absence.Rtab \\
        --annotations ${annotations} \\
        --cog-sets ${cog_sets}
    """

    stub:
    """
    mkdir distances
    """
}

process PICK_GENOMES {
    input:
    path distances
    path config          // picking YAML
    path phylogroups     // optional: [] if not given
    path extra_genomes   // optional: [] if not given

    output:
    path 'picks.tsv', emit: picks

    script:
    def pg = phylogroups ? "--phylogroups ${phylogroups}" : ''
    def extra = extra_genomes ? "--extra-genomes ${extra_genomes}" : ''
    """
    pick_genomes.py --distances ${distances} --config ${config} ${pg} ${extra}
    """

    stub:
    """
    printf 'set\\tlabel\\tgenome_id\\n' > picks.tsv
    """
}

process PLOT_PCOA {
    input:
    path distances
    path picks
    path phylogroups     // optional: [] if not given
    path extra_genomes   // optional: [] if not given

    output:
    path 'plots', emit: plots
    path 'plots_supp', emit: plots_supp   // variance explained per component, for the supplementary

    script:
    def pg = phylogroups ? "--phylogroups ${phylogroups}" : ''
    def extra = extra_genomes ? "--extra-genomes ${extra_genomes}" : ''
    """
    plot_pcoa.py --distances ${distances} --picks ${picks} ${pg} ${extra}
    """

    stub:
    """
    mkdir plots plots_supp
    """
}

process SUMMARISE_PICKS {
    input:
    path picks
    path pangenome      // Panaroo output directory
    path qc_report      // qc_report.tsv (genome lengths)
    path distances

    output:
    path 'picks_summary.tsv', emit: summary

    script:
    """
    summarise_picks.py \\
        --picks ${picks} \\
        --presence ${pangenome}/gene_presence_absence.Rtab \\
        --qc-report ${qc_report} \\
        --distances ${distances}
    """

    stub:
    """
    touch picks_summary.tsv
    """
}
