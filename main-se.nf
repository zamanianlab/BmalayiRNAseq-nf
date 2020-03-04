#!/usr/bin/env nextflow

// Edit nextflow.configuration!
data=config.brc_location
aux=config.aux_location
output=config.output_location

large_core=config.large_core
small_core=config.small_core


params.dir = "181220_BCCV00ANXX"

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
       set id, file("${id}_trim.fq.gz") into trimmed_reads
      // set val(id_out), file(id_out) into fq_trim
       file("*_trimout.txt") into trim_log


   script:
   id_out = id.replace('.fastq.gz', '_trim.fq.gz')


   """
       trimmomatic SE -threads ${large_core} ${reads} -baseout ${id}.fq.gz ILLUMINACLIP:/home/linuxbrew/.linuxbrew/Cellar/trimmomatic/0.36/share/trimmomatic/adapters/TruSeq3-PE.fa:2:80:10 MINLEN:50 &> ${id}_trimout.txt
       trimmomatic SE -phred33 -threads ${large_core} ${reads} ${name_out} ILLUMINACLIP:${adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15 &> ${reads}_trimout.txt
       rm ${id}_1U.fq.gz
       rm ${id}_2U.fq.gz
   """
}
