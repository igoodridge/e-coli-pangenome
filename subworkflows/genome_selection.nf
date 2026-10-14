include { GENOME_DISTANCES; PICK_GENOMES; PLOT_PCOA; SUMMARISE_PICKS } from '../modules/05_genome_selection'

// Jaccard distances and PCoA per COG set, the picked genomes, plots, and a summary of the picks.
// Distances only rerun when the pangenome, annotations or COG sets change; picks and plots when picking.yml does.
workflow GENOME_SELECTION {
    take:
    pangenome     // path: Panaroo output directory
    annotations   // path: gene_annotations.tsv
    qc_report     // path: qc_report.tsv (genome lengths)
    opts          // map: pipeline params (cog_sets, picking, phylogroups, extra_genomes)

    main:
    phylogroups = opts.phylogroups ?: []
    extra_genomes = opts.extra_genomes ?: []

    GENOME_DISTANCES(pangenome, annotations, opts.cog_sets)
    PICK_GENOMES(GENOME_DISTANCES.out.distances, opts.picking, phylogroups, extra_genomes)
    PLOT_PCOA(GENOME_DISTANCES.out.distances, PICK_GENOMES.out.picks, phylogroups, extra_genomes)
    SUMMARISE_PICKS(PICK_GENOMES.out.picks, pangenome, qc_report, GENOME_DISTANCES.out.distances)

    emit:
    distances = GENOME_DISTANCES.out.distances
    picks     = PICK_GENOMES.out.picks
    plots     = PLOT_PCOA.out.plots.mix(PLOT_PCOA.out.plots_supp)
    summary   = SUMMARISE_PICKS.out.summary
}
