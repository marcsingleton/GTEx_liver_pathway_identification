#!/usr/bin/env nextflow

// IO paths
params.output_path = "$projectDir/results/"
params.scripts_path = "$projectDir/code/scripts/"

// Data paths
params.reference_gtf = "$projectDir/data/GTEx/gene_model/references_v8_gencode.v26.GRCh38.genes.gtf"
params.GTEx_bulk_count = "$projectDir/data/GTEx/bulk_count/"

process protein_coding{
    publishDir "$params.output_path/protein_coding/"

    input:
        path gtf_path
    
    output:
        path "genes.tsv"

    script:
    """
    python $params.scripts_path/protein_coding.py $gtf_path
    """
}

process exploratory{
    publishDir "$params.output_path/exploratory/"

    input:
        path gtf_path
        path read_path
    
    output:
        path "*.png"
    
    script:
    """
    python $params.scripts_path/exploratory.py $gtf_path $read_path
    """
}

workflow{
    reference_gtf = Channel.fromPath(params.reference_gtf)
    protein_coding_gtf = protein_coding(reference_gtf)
    exploratory(protein_coding_gtf, Channel.fromPath(params.GTEx_bulk_count))
}
