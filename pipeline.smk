# IO paths
output_path = 'results_smk/'
scripts_path = 'code/scripts/'

# Data paths
reference_gtf = 'data/GTEx/gene_model/references_v8_gencode.v26.GRCh38.genes.gtf'
GTEx_bulk_count = 'data/GTEx/bulk_count/'

# Extract protein coding genes and perform initial exploratory analysis
rule protein_coding:
    input:
        gtf_path = reference_gtf
    output:
        protein_coding_gtf = f'{output_path}/protein_coding/genes.tsv'
    shell:
        f'''
        python {scripts_path}/protein_coding.py {{input.gtf_path}} -o {output_path}/protein_coding/
        '''

rule exploratory:
    input:
        gtf_path = rules.protein_coding.output.protein_coding_gtf,
        read_path = GTEx_bulk_count
    shell:
        f'''
        python {scripts_path}/exploratory.py {{input.gtf_path}} {{input.read_path}} -o {output_path}/exploratory/
        '''

# Compute and plot differential gene analysis
rule differential_compute:
    input:
        gtf_path = rules.protein_coding.output.protein_coding_gtf,
        read_path = GTEx_bulk_count
    output:
        anova = f'{output_path}/differential/unpooled/anova.tsv',
        ttest = f'{output_path}/differential/pooled/ttest.tsv',
        corr = f'{output_path}/differential/pooled/corr.tsv'
    shell:
        f'''
        python {scripts_path}/differential_compute.py {{input.gtf_path}} {{input.read_path}} -o {output_path}/differential/
        '''

rule differential_plot:
    input:
        anova = rules.differential_compute.output.anova,
        ttest = rules.differential_compute.output.ttest,
        corr = rules.differential_compute.output.corr,
    shell:
        f'''
        python {scripts_path}/differential_plot.py {{input.anova}} {{input.ttest}} {{input.corr}} -o {output_path}/differential/
        '''

# Compute and plot network analysis
rule network_compute:
    input:
        corr = rules.differential_compute.output.corr
    output:
        graph = f'{output_path}/network/graph.tsv',
        components = f'{output_path}/network/components.tsv'
    shell:
        f'''
        python {scripts_path}/network_compute.py {{input.corr}} -o {output_path}/network/
        '''

rule network_plot:
    input:
        corr = rules.differential_compute.output.corr,
        graph = rules.network_compute.output.graph,
        components = rules.network_compute.output.components
    shell:
        f'''
        python {scripts_path}/network_plot.py {{input.corr}} {{input.graph}} {{input.components}} -o {output_path}/network/
        '''
