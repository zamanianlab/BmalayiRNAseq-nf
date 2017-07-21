#!/usr/bin/env nextflow

// Edit nextflow.configuration!

aux=config.aux_location
data=config.data_location
output=config.output_location

large_core=config.large_core
small_core=config.small_core

// ** - Get txt file of SRA accession IDs from 'auxillary' folder
sra_file = Channel.fromPath(aux + "SRR_Acc_List.txt")

// ** - Download SRA files based on text file list of SRA accession IDs (goes to ncbi folder)
process fetch_SRA {
    
    input:
        file("SRR_Acc_List.txt") from sra_file

    script:

    sra_list="SRR_Acc_List.txt"

    """ 

    while read line     
    do           
        echo \$line
        prefetch \$line 
    done <${sra_list} 

    """
}