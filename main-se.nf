#!/usr/bin/env nextflow

// Edit nextflow.configuration!
data=config.brc_location
aux=config.aux_location
output=config.output_location

large_core=config.large_core
small_core=config.small_core

// Additional params (call from commandline: --dir "181220_BCCV00ANXX")
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
       set id_out, file("${id_out}") into trimmed_fqs
      // set val(id_out), file(id_out) into fq_trim
       file("*_trimout.txt") into trim_log


   script:
   id_out = id.replace('.fastq.gz', '_trim.fq.gz')


   """
       trimmomatic SE -threads ${large_core} ${reads} ${id_out} ILLUMINACLIP:/home/linuxbrew/.linuxbrew/Cellar/trimmomatic/0.36/share/trimmomatic/adapters/TruSeq3-PE.fa:2:80:10 MINLEN:50 &> ${id_out}_trimout.txt

   """
}
trimmed_fqs.into { trimmed_reads_hisat }

////////////////////////////////////////////////
// ** - Fetch Parasite (P) reference genome (fa.gz) and gene annotation file (gtf.gz)
////////////////////////////////////////////////

// release="WBPS13"
// species="brugia_malayi"
// prjn="PRJNA10729"
// prefix="ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/${release}/species/${species}/${prjn}"
//
// process fetch_reference {
//
//     publishDir "${output}/reference/", mode: 'copy'
//
//     output:
//         file("geneset.gtf.gz") into geneset_gtf
//         file("reference.fa.gz") into reference_fa
//
//     """
//         echo '${prefix}'
//         curl ${prefix}/${species}.${prjn}.${release}.canonical_geneset.gtf.gz > geneset.gtf.gz
//         curl ${prefix}/${species}.${prjn}.${release}.genomic.fa.gz > reference.fa.gz
//
//     """
// }
// geneset_gtf.into { geneset_hisat; geneset_stringtie }
// reference_fa.into { reference_hisat; reference_bwa}
