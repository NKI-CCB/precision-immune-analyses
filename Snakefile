rule all:
    input:
        'plots/boxplot_stroma_all_density_clinical_variables.pdf',
        'plots/boxplot_dcis_all_density_clinical_variables.pdf',
        'plots/boxplot_stroma_fibrosis_density_clinical_variables.pdf',
        'plots/boxplot_dcis_fibrosis_density_clinical_variables.pdf',
        'plots/boxplot_stroma_erpos_her2neg_density_clinical_variables.pdf',
        'plots/boxplot_dcis_erpos_her2neg_density_clinical_variables.pdf',

rule boxplot_dens_clin_stroma:
    output:
        plot='plots/boxplot_stroma_{subset}_density_clinical_variables.pdf',
        tests='results/significance_stroma_{subset}.xlsx',
        stats='results/summary_stats_stroma_{subset}.xlsx',
    input:
        density='data/cell_density_Stroma.tsv',
        density_ki67='data/ki67_density.tsv',
        clinical='data/clin_DBL_v_str_dens.tsv',
        script='src/beeswarm_clin_dens.R',
    shell:
        'Rscript {input.script} {input.density} {input.density_ki67} Tissue {input.clinical}'
        ' {wildcards.subset} {output.plot} {output.stats} {output.tests}'

rule boxplot_dens_clin_dcis:
    output:
        plot='plots/boxplot_dcis_{subset}_density_clinical_variables.pdf',
        tests='results/significance_dcis_{subset}.xlsx',
        stats='results/summary_stats_dcis_{subset}.xlsx',
    input:
        density='data/cell_density_DCIS.tsv',
        density_ki67='data/ki67_density.tsv',
        clinical='data/clin_DBL_v_str_dens.tsv',
        script='src/beeswarm_clin_dens.R',
    shell:
        'Rscript {input.script} {input.density} {input.density_ki67} DCIS {input.clinical}'
        ' {wildcards.subset} {output.plot} {output.stats} {output.tests}'
