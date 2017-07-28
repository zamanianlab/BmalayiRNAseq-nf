#!/usr/bin/env nextflow

// Edit nextflow.configuration!

aux=config.aux_location
data=config.data_location
output=config.output_location

large_core=config.large_core
small_core=config.small_core

// ** - Get txt file of SRA accession IDs from 'auxillary' folder
sra_file = Channel.fromPath(aux + "SRR_Acc_List.txt")


// ** - Recurse through subdirectories to get all fastqs
// fq_set = Channel.fromPath(data + "fq/*.fastq.gz")
//                 .map { n -> [ n.getName(), n ] }

// Fetch fqs; alternative suffixes
Channel.fromFilePairs(data +'fq/*_{1,2}.fastq.gz', flat: true)
        .into { read_pairs }

// SKIP TRIMMING (READS ARE ALREADY TRIMMED)
// process trim {

//     tag { fq_id }

//     publishDir "${data}/fq_trim/", mode: 'move'

//     input:
//         set fq_id, file(forward), file(reverse) from read_pairs

//     output:
//         set file("${fq_id}_1P.fq.gz"), file("${fq_id}_2P.fq.gz") into trim_output
//         //file "${fq_id}.trim_log.txt" into trim_logs

//     """
//     trimmomatic PE -threads ${large_core} $forward $reverse -baseout ${fq_id}.fq.gz ILLUMINACLIP:/home/linuxbrew/.linuxbrew/Cellar/trimmomatic/0.36/share/trimmomatic/adapters/TruSeq3-PE.fa:2:80:10 MINLEN:75 
//     rm ${fq_id}_1U.fq.gz
//     rm ${fq_id}_2U.fq.gz
//     """

// }


// process trimmomatic {

//     cpus small_core

//     tag { name }

//     input:
//         set val(name), file(reads) from fq_set

//     output:
//         file(name_out) into trimmed_reads

//     script:
//     name_out = name.replace('.fastq.gz', '_trim.fq.gz')

//     """
//         trimmomatic SE -phred33 -threads ${small_core} ${reads} ${name_out} ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15 
//     """
// }


// ** - Fetch reference genome (fa.gz) and gene annotation file (gtf.gz)
release="WBPS9"
species="brugia_malayi"
prjn="PRJNA10729"
prefix="ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/${release}/species/${species}/${prjn}"

process fetch_reference {

    publishDir "${data}/reference/", mode: 'copy'
    
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
extract_exons_py = file("${aux}/scripts/hisat2_extract_exons.py")
extract_splice_py = file("${aux}/scripts/hisat2_extract_splice_sites.py")

process hisat2_indexing {

   publishDir "${data}/reference/", mode: 'copy'

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

    publishDir "${data}/reference/", mode: 'copy'

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


// ** - ALIGNMENT
process align {

    cpus small_core

    tag { srid }

    input:
        set val(srid), file(forward), file(reverse) from read_pairs
        file hs2_indices from hs2_indices.first()

    output:
        set val(srid), file("${srid}.bam"), file("${srid}.bam.bai") into hisat2_bams
        file "${prefix}.hisat2_log.txt" into alignment_logs

    script:
        index_base = hs2_indices[0].toString() - ~/.\d.ht2/

    """
        hisat2 -p ${small_core} -x $index_base -1 ${forward} -2 ${reverse} -S ${srid}.sam --rg-id "${srid}" --rg "SM:${srid}" --rg "PL:ILLUMINA" 2> ${srid}.hisat2_log.txt
        samtools view -bS ${srid}.sam > ${srid}.unsorted.bam
        samtools flagstat ${srid}.unsorted.bam
        samtools sort -@ ${small_core} -o ${srid}.bam ${srid}.unsorted.bam
        samtools index -b ${srid}.bam
        rm *sam
        rm *unsorted.bam

    """
}



process stringtie_counts {

    publishDir "output/expression", mode: 'copy'

    cpus small_core

    tag { srid }

    input:
        set val(srid), file(bam), file(bai) from hisat2_bams
        file("geneset.gtf.gz") from geneset_stringtie.first()

    output:
        file("${srid}/*") into stringtie_exp

    """ 
        zcat geneset.gtf.gz > geneset.gtf
        stringtie -p ${small_core} -G geneset.gtf -A ${srid}/${srid}_abund.tab -e -B -o ${srid}/${srid}_expressed.gtf ${bam}
    """
}



// prepDE = file("auxillary/scripts/prepDE.py")

// process stringtie_table_counts {

//     echo true

//     publishDir "output/diffexp", mode: 'copy'

//     cpus small_core

//     tag { sample_id }

//     input:
//         val(sample_file) from stringtie_exp.toSortedList()

//     output:
//         file ("gene_count_matrix.csv") into gene_count_matrix
//         file ("transcript_count_matrix.csv") into transcript_count_matrix

//     """
//         for i in ${sample_file.flatten().join(" ")}; do
//             bn=`basename \${i}`
//             full_path=`dirname \${i}`
//             sample_name=\${full_path##*/}
//             echo "\${sample_name} \${i}"
//             mkdir -p expression/\${sample_name}
//             ln -s \${i} expression/\${sample_name}/\${bn}
//         done;
//         python ${prepDE} -i expression -l 50 -g gene_count_matrix.csv -t transcript_count_matrix.csv

//     """
// }
