#!/usr/bin/env nextflow

// Edit nextflow.configuration!
aux=config.aux_location
data=config.data_location
output=config.output_location

large_core=config.large_core
small_core=config.small_core

// ** - Get txt file of SRA accession IDs from 'auxillary' folder
//sra_file = Channel.fromPath(aux + "SRR_Acc_List_UGA.txt")

////////////////////////////////////////////////
// ** - Pull in fq files (paired)
////////////////////////////////////////////////

Channel.fromFilePairs(data +'fq/*_{1,2}.fastq.gz', flat: true)
        .into { read_pairs }

////////////////////////////////////////////////
// ** TRIM READS
////////////////////////////////////////////////

process trimmomatic {

   cpus large_core
   tag { id }
   publishDir "${output}/trim_stats/", mode: 'copy', pattern: '*_trimout.txt'

   input:
       set val(id), file(forward), file(reverse) from read_pairs

   output:
       set id, file("${id}_1P.fq.gz"), file("${id}_2P.fq.gz") into trimmed_read_pairs
       file("*_trimout.txt") into trim_log

   script:
   id_out = id.replace('.fastq.gz', '_trim.fq.gz')

   """
       trimmomatic PE -threads ${large_core} $forward $reverse -baseout ${id}.fq.gz ILLUMINACLIP:/home/linuxbrew/.linuxbrew/Cellar/trimmomatic/0.36/share/trimmomatic/adapters/TruSeq3-PE.fa:2:80:10 MINLEN:50 &> ${id}_trimout.txt
       rm ${id}_1U.fq.gz
       rm ${id}_2U.fq.gz
   """
}


////////////////////////////////////////////////
// ** - Fetch Parasite (P) reference genome (fa.gz) and gene annotation file (gtf.gz)
////////////////////////////////////////////////

release="WBPS13"
species="brugia_malayi"
prjn="PRJNA10729"
prefix="ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/${release}/species/${species}/${prjn}"

process fetch_reference {

    publishDir "${output}/reference/", mode: 'copy'

    output:
        file("geneset.gtf.gz") into geneset_gtf
        file("reference.fa.gz") into reference_hisat

    """
        echo '${prefix}'
        curl ${prefix}/${species}.${prjn}.${release}.canonical_geneset.gtf.gz > geneset.gtf.gz
        curl ${prefix}/${species}.${prjn}.${release}.genomic.fa.gz > reference.fa.gz

    """

}
geneset_gtf.into { geneset_hisat; geneset_stringtie }

// ** - Create HiSat2 Index using reference genome and annotation file
extract_exons = file("${aux}/scripts/hisat2_extract_exons.py")
extract_splice = file("${aux}/scripts/hisat2_extract_splice_sites.py")

process hisat2_indexing {

   publishDir "${output}/reference/", mode: 'copy'

    input:
        file("geneset.gtf.gz") from geneset_hisat
        file("reference.fa.gz") from reference_hisat

    output:
        file("splice.ss") into splice_hisat
        file("exon.exon") into exon_hisat
        file("reference.fa.gz") into reference_build_hisat

    """
        zcat geneset.gtf.gz | python ${extract_splice} - > splice.ss
        zcat geneset.gtf.gz | python ${extract_exons} - > exon.exon
    """

}

process build_hisat_index {

    publishDir "${output}/reference/", mode: 'copy'

    cpus large_core

    input:
        file("splice.ss") from splice_hisat
        file("exon.exon") from exon_hisat
        file("reference.fa.gz") from reference_build_hisat

    output:
        file "*.ht2" into hs2_indices

    """
        zcat reference.fa.gz > reference.fa
        hisat2-build -p ${large_core} --ss splice.ss --exon exon.exon reference.fa reference.hisat2_index
    """

}


////////////////////////////////////////////////
// ** - ALIGNMENT AND STRINGTIE (combined) against Parasite
////////////////////////////////////////////////

process hisat2_stringtie {

    publishDir "${output}/expression", mode: 'copy'

    cpus large_core
    tag { id }

    input:
        set val(id), file(forward), file(reverse) from trimmed_read_pairs
        file("geneset.gtf.gz") from geneset_stringtie
        file hs2_indices from hs2_indices.first()

    output:
        file "${id}.hisat2_log.txt" into alignment_logs
        file("${id}/*") into stringtie_exp

    script:
        index_base = hs2_indices[0].toString() - ~/.\d.ht2/

    """
        hisat2 -p ${large_core} -x $index_base -1 ${forward} -2 ${reverse} -S ${id}.sam --rg-id "${id}" --rg "SM:${id}" --rg "PL:ILLUMINA" 2> ${id}.hisat2_log.txt
        samtools view -bS ${id}.sam > ${id}.unsorted.bam
        rm *.sam
        samtools flagstat ${id}.unsorted.bam
        samtools sort -@ ${large_core} -o ${id}.bam ${id}.unsorted.bam
        rm *.unsorted.bam
        samtools index -b ${id}.bam
        zcat geneset.gtf.gz > geneset.gtf
        stringtie ${id}.bam -p ${large_core} -G geneset.gtf -A ${id}/${id}_abund.tab -e -B -o ${id}/${id}_expressed.gtf
        rm *.bam
        rm *.bam.bai
        rm *.gtf
    """
}



//run last
////////////////////////////////////////////////
// ** - STRINGTIE table counts
////////////////////////////////////////////////

// prepDE = file("${aux}/scripts/prepDE.py")
// process stringtie_table_counts {
//
//     echo true
//
//     publishDir "${output}/diffexp", mode: 'copy'
//
//     cpus small_core
//
//     output:
//         file ("gene_count_matrix.csv") into gene_count_matrix
//         file ("transcript_count_matrix.csv") into transcript_count_matrix
//
//     """
//         python ${prepDE} -i ${output}/expression -l 100 -g gene_count_matrix.csv -t transcript_count_matrix.csv
//
//     """
// }
