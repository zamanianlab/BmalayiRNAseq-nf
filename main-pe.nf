#!/usr/bin/env nextflow

// Edit nextflow.configuration!
data=config.brc_location
aux=config.aux_location
output=config.output_location

large_core=config.large_core
small_core=config.small_core

// Additional params (--dir "200217_AHNHN3DMXX")
params.dir = "200311_AHNNC7DMXX"
// flag for final stringtie_table_counts process (--stc)
params.stc = false

//WB genome information
params.release = "WBPS13"
params.species = "dirofilaria_immitis" //brugia_malayi
params.prjn = "PRJEB1797" //PRJNA10729

////////////////////////////////////////////////
// ** - Pull in fq files (paired)
////////////////////////////////////////////////

Channel.fromFilePairs(data + "${params.dir}/*_R{1,2}_001.fastq.gz", flat: true)
        .set { fq_pairs }

////////////////////////////////////////////////
// ** TRIM READS
////////////////////////////////////////////////

process trimmomatic {

   cpus large_core
   tag { id }
   publishDir "${output}/trim_stats/", mode: 'copy', pattern: '*_trimout.txt'

   input:
       set val(id), file(forward), file(reverse) from fq_pairs

   output:
       set id, file("${id}_1P.fq.gz"), file("${id}_2P.fq.gz") into trimmed_fq_pairs
       file("*_trimout.txt") into trim_log

   """
       trimmomatic PE -threads ${large_core} $forward $reverse -baseout ${id}.fq.gz ILLUMINACLIP:/home/linuxbrew/.linuxbrew/Cellar/trimmomatic/0.36/share/trimmomatic/adapters/TruSeq3-PE.fa:2:80:10 MINLEN:50 &> ${id}_trimout.txt
       rm ${id}_1U.fq.gz
       rm ${id}_2U.fq.gz
   """
}
trimmed_fq_pairs.set { trimmed_reads_hisat }

////////////////////////////////////////////////
// ** - Fetch genome (fa.gz) and gene annotation file (gtf.gz)
////////////////////////////////////////////////

process fetch_genome {

    publishDir "${output}/reference/", mode: 'copy'

    output:
        file("geneset.gtf.gz") into geneset_gtf
        file("reference.fa.gz") into reference_fa

    script:

        prefix="ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/${params.release}/species/${params.species}/${params.prjn}"

    """
        echo '${prefix}'
        curl ${prefix}/${params.species}.${params.prjn}.${params.release}.canonical_geneset.gtf.gz > geneset.gtf.gz
        curl ${prefix}/${params.species}.${params.prjn}.${params.release}.genomic.fa.gz > reference.fa.gz
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

    publishDir "${output}/expression", mode: 'copy', pattern: '**/*'
    publishDir "${output}/expression", mode: 'copy', pattern: '*.hisat2_log.txt'
    publishDir "${output}/bams", mode: 'copy', pattern: '*.bam'
    publishDir "${output}/bams", mode: 'copy', pattern: '*.bam.bai'

    cpus large_core
    tag { id }

    input:
        set val(id), file(forward), file(reverse) from trimmed_reads_hisat
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
        hisat2 -p ${large_core} -x $index_base -1 ${forward} -2 ${reverse} -S ${id}.sam --rg-id "${id}" --rg "SM:${id}" --rg "PL:ILLUMINA" 2> ${id}.hisat2_log.txt
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

////////////////////////////////////////////////
// ** - STRINGTIE table counts (run last with --stc flag)
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
