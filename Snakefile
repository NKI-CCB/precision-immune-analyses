wildcard_constraints:
    subset='all|erpos_her2neg|fibrosis'

rule all:
    input:
        'plots/boxplot_stroma_all_density_clinical_variables.pdf',
        'plots/boxplot_dcis_all_density_clinical_variables.pdf',
        'plots/boxplot_stroma_fibrosis_density_clinical_variables.pdf',
        'plots/boxplot_dcis_fibrosis_density_clinical_variables.pdf',
        'plots/boxplot_stroma_erpos_her2neg_density_clinical_variables.pdf',
        'plots/boxplot_dcis_erpos_her2neg_density_clinical_variables.pdf',
        'results/significance_stroma_ratio_all.xlsx',
        'results/significance_dcis_ratio_all.xlsx',
        'plots/heatmap_density_Stroma.pdf',
        'plots/heatmap_density_DCIS.pdf',
        'results/significance_til_density_all.xlsx',

rule boxplot_dens_clin_stroma:
    output:
        plot='plots/boxplot_stroma_{subset}_density_clinical_variables.pdf',
        tests='results/significance_stroma_{subset}.xlsx',
        stats='results/summary_stats_stroma_{subset}.xlsx',
    input:
        density='data/cell_density_Stroma.tsv',
        clinical='data/clin_DBL_v_str_dens.tsv',
        script='src/beeswarm_clin_dens.R',
    shell:
        'Rscript {input.script} {input.density} Tissue {input.clinical}'
        ' {wildcards.subset} {output.plot} {output.stats} {output.tests}'

rule boxplot_dens_clin_dcis:
    output:
        plot='plots/boxplot_dcis_{subset}_density_clinical_variables.pdf',
        tests='results/significance_dcis_{subset}.xlsx',
        stats='results/summary_stats_dcis_{subset}.xlsx',
    input:
        density='data/cell_density_DCIS.tsv',
        clinical='data/clin_DBL_v_str_dens.tsv',
        script='src/beeswarm_clin_dens.R',
    shell:
        'Rscript {input.script} {input.density} DCIS {input.clinical}'
        ' {wildcards.subset} {output.plot} {output.stats} {output.tests}'

rule ratios_clin_stroma:
    output:
        tests='results/significance_stroma_ratio_{subset}.xlsx',
        stats='results/summary_stats_stroma_ratio_{subset}.xlsx',
    input:
        density='data/cell_density_Stroma.tsv',
        clinical='data/clin_DBL_v_str_dens.tsv',
        script='src/ratios.R',
    shell:
        'Rscript {input.script} {input.density} Tissue {input.clinical}'
        ' {wildcards.subset} {output.stats} {output.tests}'

rule ratios_clin_dcis:
    output:
        tests='results/significance_dcis_ratio_{subset}.xlsx',
        stats='results/summary_stats_dcis_ratio_{subset}.xlsx',
    input:
        density='data/cell_density_DCIS.tsv',
        clinical='data/clin_DBL_v_str_dens.tsv',
        script='src/ratios.R',
    shell:
        'Rscript {input.script} {input.density} DCIS {input.clinical}'
        ' {wildcards.subset} {output.stats} {output.tests}'

rule dens_grade_test:
    output:
        tests='results/significance_stroma_densgrade_{subset}.xlsx',
    input:
        density='data/cell_density_Stroma.tsv',
        clinical='data/clin_DBL_v_str_dens.tsv',
        script='src/dens_and_grade.R',
    shell:
        'Rscript {input.script} {input.density} Tissue {input.clinical}'
        ' {wildcards.subset} {output.tests}'


rule pathologist_scored_tils:
    output:
        tests='results/significance_til_density_{subset}.xlsx',
    input:
        density='data/cell_density_Stroma.tsv',
        clinical='data/clin_DBL_v_str_dens.tsv',
        til_scores='data/20220324_TILscores_revised.tsv',
        script='src/til_stats.R',
    shell:
        'Rscript {input.script} {input.density} Tissue {input.clinical}'
        '  {wildcards.subset} {input.til_scores} {output.tests}'

rule clustering:
    output:
       'plots/heatmap_density_Stroma.pdf',
       'plots/heatmap_density_DCIS.pdf',
    input:
        density_stroma='data/cell_density_Stroma.tsv',
        density_dcis='data/cell_density_DCIS.tsv',
        clinical='data/clin_DBL_v_str_dens.tsv',
        script='src/clustering.R',
    shell:
        'Rscript {input.script}'
