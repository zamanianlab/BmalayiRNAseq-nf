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


// ** - Covert SRA files to fastqs
// process fetch_reads {

//     publishDir "${data}/fq/", mode: 'copy'
    
//     input:
//         file("SRR_Acc_List.txt") from sra_file

//     output:
//         file("*")

//     script:

//     sra_list="SRR_Acc_List.txt"

//     """ 

//     while read line     
//     do           
//         touch \$line.txt
//     done <${sra_list} 

//     """
// //**  fastq-dump --gzip \$line  
// }


// ** - Recurse through subdirectories to get all fastqs
// fq_set = Channel.fromPath(data + "sra/*.fastq.gz")
//                 .map { n -> [ n.getName(), n ] }


** - Fetch reference genome (fa.gz) and gene annotation file (gtf.gz)
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
        echo '${release}'
        echo '${species}'
        echo '${prefix}'
        wget -nc -r -nH --cut-dirs=7 --no-parent --reject="index.html*" -A 'canonical_geneset.gtf.gz','genomic.fa.gz' $prefix
        mv '${species}/${prjn}/${species}.${prjn}.${release}.canonical_geneset.gtf.gz' geneset.gtf.gz
        mv '${species}/${prjn}/${species}.${prjn}.${release}.genomic.fa.gz' reference.fa.gz

    """
}
geneset_gtf.into { geneset_hisat; geneset_stringtie }


// ** - Create HiSat2 Index using reference genome and annotation file

// extract_exons_py = file("${aux}/scripts/hisat2_extract_exons.py")
// extract_splice_py = file("${aux}/scripts/hisat2_extract_splice_sites.py")

// process hisat2_indexing {

//     input:
//         file("geneset.gtf.gz") from geneset_hisat
//         file("reference.fa.gz") from reference_hisat

//     output:
//         file("splice.ss") into splice_hisat
//         file("exon.exon") into exon_hisat
//         file("reference.fa.gz") into reference_build_hisat

//     """
//         zcat geneset.gtf.gz | python ${extract_splice_py} - > splice.ss
//         zcat geneset.gtf.gz | python ${extract_exons_py} - > exon.exon
//     """

// }

// process build_hisat_index {

//     cpus small_core

//     input:
//         file("splice.ss") from splice_hisat
//         file("exon.exon") from exon_hisat
//         file("reference.fa.gz") from reference_build_hisat

//     output:
//         file "*.ht2" into hs2_indices

//     """
//         zcat reference.fa.gz > reference.fa
//         hisat2-build -p ${small_core} --ss splice.ss --exon exon.exon reference.fa reference.hisat2_index
//     """

// }

// process align {

//     cpus large_core

//     publishDir "${output}/bam", mode: 'copy'

//     tag { reads }

//     input:
//         set val(name), file(reads) from fq_set
//         file hs2_indices from hs2_indices.first()

//     output:
//         set val("${prefix}"), file("${prefix}.bam"), file("${prefix}.bam.bai") into hisat2_bams
//         file "${prefix}.hisat2_log.txt" into alignment_logs

//     script:
//         index_base = hs2_indices[0].toString() - ~/.\d.ht2/
//         prefix = reads[0].toString() - ~/(\.fastq\.gz)$/

//     """ 
//         hisat2 -p ${large_core} -x $index_base -U ${reads} -S ${prefix}.sam 2> ${prefix}.hisat2_log.txt
//         samtools view -bS ${prefix}.sam > ${prefix}.unsorted.bam
//         samtools flagstat ${prefix}.unsorted.bam
//         samtools sort -@ ${small_core} -o ${prefix}.bam ${prefix}.unsorted.bam
//         samtools index -b ${prefix}.bam
//         rm *sam
//         rm *unsorted.bam
//     """
// }


// process stringtie_counts {

//     publishDir "output/expression", mode: 'copy'

//     cpus small_core

//     tag { srid }

//     input:
//         set val(srid), file(bam), file(bai) from hisat2_bams
//         file("geneset.gtf.gz") from geneset_stringtie.first()

//     output:
//         file("${srid}/*") into stringtie_exp

//     """ 
//         zcat geneset.gtf.gz > geneset.gtf
//         stringtie -p ${small_core} -G geneset.gtf -A ${srid}/${srid}_abund.tab -e -B -o ${srid}/${srid}_expressed.gtf ${bam}
//     """
// }



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
