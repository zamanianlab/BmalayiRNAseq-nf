#!/usr/bin/env nextflow

// Edit nextflow.configuration!

data_location=config.data_location
aux_location=config.aux_location
large_core=config.large_core
small_core=config.small_core

// ** - Get txt file of SRA accession IDs from 'auxillary' folder
sra_file = Channel.fromPath(aux_location + "SRR_Acc_List.txt")

// ** - Download SRA files based on text file list of SRA accession IDs
process fetch_reads {

    publishDir "data/SRA/", mode: 'copy', pattern: 'meta.pipeline.txt'
    
    input:
        file("SRR_Acc_List.txt") from sra_file

    output:
        file("meta.pipeline.txt")
        file("*.fastq.gz")

    script:

    sra_list="SRR_Acc_List.txt"

    """ 

    echo ${sra_list}
    while read line     
    do           
        echo \$line >> meta.pipeline.txt
        fastq-dump --gzip \$line 
    done <${sra_list} 

    """
}

// ** - Recurse through subdirectories to get all fastqs
fq_set = Channel.fromPath(data_location + "**/*.fastq.gz")
                .map { n -> [ n.getName(), n ] }