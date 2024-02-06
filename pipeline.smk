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
