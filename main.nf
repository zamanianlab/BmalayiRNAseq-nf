#!/usr/bin/env nextflow

// Edit nextflow.configuration!

data_location=config.data_location
aux_location=config.aux_location
large_core=config.large_core
small_core=config.small_core

// ** - Get txt file of SRA accession IDs from 'auxillary' folder
sra_file = Channel.fromPath(aux_location + "SRR_Acc_List_sample.txt")

// ** - Download SRA files based on text file list of SRA accession IDs
process fetch_reads {

    publishDir "data/SRA/", mode: 'copy', pattern: 'meta.pipeline.txt'
    
    input:
        file("SRR_Acc_List_sample.txt") from sra_file

    output:
        file("meta.pipeline.txt")
        file("*.fastq.gz")

    script:

    sra_list="SRR_Acc_List_sample.txt"

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



// ** - Fetch reference genome (fa.gz) and gene annotation file (gtf.gz)
release="WBPS9"
species="brugia_malayi"
prjn="PRJNA10729"
prefix="ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/${release}/species/${species}/${prjn}"

process fetch_reference {

    publishDir "data/reference/", mode: 'copy'
    
    output:
        file("${species}/${prjn}/${species}.${prjn}.${release}.canonical_geneset.gtf.gz") into geneset_gtf
        file("${species}/${prjn}/${species}.${prjn}.${release}.genomic.fa.gz") into reference_hisat

    """
        echo '${release}'
        echo '${species}'
        echo '${prefix}'
        wget -nc -r -nH --cut-dirs=7 --no-parent --reject="index.html*" -A 'canonical_geneset.gtf.gz','genomic.fa.gz' $prefix

    """

}
geneset_gtf.into { geneset_hisat; geneset_stringtie }


// ** - Create HiSat2 Index using reference genome and annotation file

extract_exons_py = file("auxillary/scripts/hisat2_extract_exons.py")
extract_splice_py = file("auxillary/scripts/hisat2_extract_splice_sites.py")

process hisat2_indexing {

    input:
        file("geneset.gtf.gz") from geneset_hisat
        file("reference.fa.gz") from reference_hisat

    output:
        file("splice.ss") into splice_hisat
        file("exon.exon") into exon_hisat
        file("reference.fa.gz") into reference_build_hisat

    """
        zcat geneset.gtf.gz | python ${extract_splice_py} - > splice.ss
        zcat geneset.gtf.gz | python ${extract_exons_py} - > exon.exon
    """

}

process build_hisat_index {

    cpus small_core

    input:
        file("splice.ss") from splice_hisat
        file("exon.exon") from exon_hisat
        file("reference.fa.gz") from reference_build_hisat

    output:
        file "*.ht2" into hs2_indices

    """
        zcat reference.fa.gz > reference.fa
        hisat2-build -p ${small_core} --ss splice.ss --exon exon.exon reference.fa reference.hisat2_index
    """

}

