rule all:
    input:
        'plots/boxplot_stroma_density_clinical_variables.pdf',
        'plots/boxplot_dcis_density_clinical_variables.pdf',

rule boxplot_dens_clin_stroma:
    output:
        plot='plots/boxplot_stroma_density_clinical_variables.pdf',
        tests='results/significance_stroma.xlsx',
        stats='results/summary_stats_stroma.xlsx',
    input:
        density='data/cell_density_Stroma.tsv',
        density_ki67='data/ki67_density.tsv',
        clinical='data/clin_DBL_v_str_dens.tsv',
        script='src/beeswarm_clin_dens.R',
    shell:
        'Rscript {input.script} {input.density} {input.density_ki67} Tissue {input.clinical}'
        ' {output.plot} {output.stats} {output.tests}'

rule boxplot_dens_clin_dcis:
    output:
        plot='plots/boxplot_dcis_density_clinical_variables.pdf',
        tests='results/significance_dcis.xlsx',
        stats='results/summary_stats_dcis.xlsx',
    input:
        density='data/cell_density_DCIS.tsv',
        density_ki67='data/ki67_density.tsv',
        clinical='data/clin_DBL_v_str_dens.tsv',
        script='src/beeswarm_clin_dens.R',
    shell:
        'Rscript {input.script} {input.density} {input.density_ki67} DCIS {input.clinical}'
        ' {output.plot} {output.stats} {output.tests}'
