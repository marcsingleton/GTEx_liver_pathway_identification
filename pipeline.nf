#!/usr/bin/env nextflow

// IO paths
params.output_path = "$projectDir/results_nf/"
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

process differential_compute{
    publishDir "$params.output_path/differential/"

    input:
        path gtf_path
        path read_path

    output:
        path "unpooled/anova.tsv", emit: anova
        path "pooled/ttest.tsv", emit: ttest
        path "pooled/corr.tsv", emit: corr

    script:
    """
    python $params.scripts_path/differential_compute.py $gtf_path $read_path
    """
}

process differential_plot{
    publishDir "$params.output_path/differential/"

    input:
        path anova
        path ttest
        path corr
    
    output:
        path "**.png"
    
    script:
    """
    python $params.scripts_path/differential_plot.py $anova $ttest $corr
    """

}

process network_compute{
    publishDir "$params.output_path/network/"

    input:
        path corr
    
    output:
        path "graph.tsv", emit: graph
        path "components.tsv", emit: components
    
    script:
    """
    python $params.scripts_path/network_compute.py $corr
    """
}

process network_plot{
    publishDir "$params.output_path/network/"

    input:
        path corr
        path graph
        path components
    
    output:
        path "**.png"
    
    script:
    """
    python $params.scripts_path/network_plot.py $corr $graph $components
    """
}

workflow{
    // Make channels from input data paths
    reference_gtf = Channel.fromPath(params.reference_gtf)
    GTEx_bulk_count = Channel.fromPath(params.GTEx_bulk_count)

    // Extract protein coding genes and perform initial exploratory analysis
    protein_coding_gtf = protein_coding(reference_gtf)
    exploratory_plots = exploratory(protein_coding_gtf, GTEx_bulk_count)
    
    // Compute and plot differential gene analysis
    diff_results = differential_compute(protein_coding_gtf, GTEx_bulk_count)
    diff_plots = differential_plot(diff_results)

    // Compute and plot network analysis
    network_results = network_compute(diff_results.corr)
    network_plots = network_plot(diff_results.corr, network_results.graph, network_results.components)
}
