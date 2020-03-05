#!/usr/bin/env nextflow

// Edit nextflow.configuration!
data=config.brc_location
aux=config.aux_location
output=config.output_location

large_core=config.large_core
small_core=config.small_core

// Additional params (call from commandline: --dir "181220_BCCV00ANXX")
params.dir = "181220_BCCV00ANXX"
params.stc = false // flag for final process (stringtie_table_counts: --stc)

////////////////////////////////////////////////
// ** - Pull in fq files (single)
////////////////////////////////////////////////

fqs = Channel.fromPath(data + "${params.dir}/*.fastq.gz")
                        .map { n -> [ n.getName(), n ] }

////////////////////////////////////////////////
// ** TRIM READS
////////////////////////////////////////////////

process trimmomatic {

   cpus large_core
   tag { id }
   publishDir "${output}/trim_stats/", mode: 'copy', pattern: '*_trimout.txt'

   input:
       set val(id), file(reads) from fqs

   output:
       set id_out, file("${id_out}_trim.fq.gz") into trimmed_fqs
      // set val(id_out), file(id_out) into fq_trim
       file("*_trimout.txt") into trim_log

   script:
   id_out = id.replace('.fastq.gz', '')

   """
       trimmomatic SE -threads ${large_core} ${reads} ${id_out}_trim.fq.gz ILLUMINACLIP:/home/linuxbrew/.linuxbrew/Cellar/trimmomatic/0.36/share/trimmomatic/adapters/TruSeq3-SE.fa:2:80:10 MINLEN:50 &> ${id_out}_trimout.txt

   """
}
trimmed_fqs.into { trimmed_reads_hisat }

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
        file("reference.fa.gz") into reference_fa

    """
        echo '${prefix}'
        curl ${prefix}/${species}.${prjn}.${release}.canonical_geneset.gtf.gz > geneset.gtf.gz
        curl ${prefix}/${species}.${prjn}.${release}.genomic.fa.gz > reference.fa.gz

    """
}
geneset_gtf.into { geneset_hisat; geneset_stringtie }
reference_fa.into { reference_hisat; reference_bwa}

////////////////////////////////////////////////
// ** - HiSat2/Stringtie pipeline
////////////////////////////////////////////////

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

// alignment and stringtie combined
process hisat2_stringtie {

    publishDir "${output}/expression", mode: 'copy', pattern: '/*'
    publishDir "${output}/expression", mode: 'copy', pattern: '*.hisat2_log.txt'
    publishDir "${output}/bams", mode: 'copy', pattern: '*.bam'
    publishDir "${output}/bams", mode: 'copy', pattern: '*.bam.bai'

    cpus large_core
    tag { id }

    input:
        set val(id), file(reads) from trimmed_reads_hisat
        file("geneset.gtf.gz") from geneset_stringtie
        file hs2_indices from hs2_indices.first()

    output:
        file "${id}.hisat2_log.txt" into alignment_logs
        file("${id}/*") into stringtie_exp
        file("${id}.bam") into bam_files
        file("${id}.bam.bai") into bam_indexes

    script:
        index_base = hs2_indices[0].toString() - ~/.\d.ht2/

    """
        hisat2 -p ${large_core} -x $index_base -U ${reads} -S ${id}.sam --rg-id "${id}" --rg "SM:${id}" --rg "PL:ILLUMINA" 2> ${id}.hisat2_log.txt
        samtools view -bS ${id}.sam > ${id}.unsorted.bam
        rm *.sam
        samtools flagstat ${id}.unsorted.bam
        samtools sort -@ ${large_core} -o ${id}.bam ${id}.unsorted.bam
        rm *.unsorted.bam
        samtools index -b ${id}.bam
        zcat geneset.gtf.gz > geneset.gtf
        stringtie ${id}.bam -p ${large_core} -G geneset.gtf -A ${id}/${id}_abund.tab -e -B -o ${id}/${id}_expressed.gtf
        rm *.gtf
    """
}

//run last with flag
////////////////////////////////////////////////
// ** - STRINGTIE table counts
////////////////////////////////////////////////

prepDE = file("${aux}/scripts/prepDE.py")
process stringtie_table_counts {

    echo true

    publishDir "${output}/counts", mode: 'copy'

    cpus small_core

    when:
      params.stc

    output:
        file ("gene_count_matrix.csv") into gene_count_matrix
        file ("transcript_count_matrix.csv") into transcript_count_matrix

    """
        python ${prepDE} -i ${output}/expression -l 150 -g gene_count_matrix.csv -t transcript_count_matrix.csv

    """
}
