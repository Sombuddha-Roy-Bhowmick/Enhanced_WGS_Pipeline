#!/usr/bin/env nextflow


// Define parameters with default values
params.threads = params.threads ?: 'default_value'
def Ref_hg38 = "Homo_sapiens_assembly38.fasta"
def germline_resource = "af-only-gnomad.hg38.vcf.gz"
def pon = "1000g_pon.hg38.vcf.gz"
def dbsnp = "Homo_sapiens_assembly38.dbsnp138.vcf"
params.workdir = "$PWD"

// Channels for normal and tumor BAM files
def normal_bam = "${params.workdir}/Normal_BQSR/normal_bqsr.bam"
def tumor_bam = "${params.workdir}/Tumor_BQSR/tumor_bqsr.bam"

def chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

// Process to run Mutect2 for each chromosome
process runMutect2 {
    input:
    val chr

    output:
    path "Mutect2/${chr}.vcf"

    script:
    """
    mkdir -p Mutect2
    source /home/sombuddha/anaconda3/etc/profile.d/conda.sh
    conda activate wgs
    gatk --java-options "-Xmx5g" Mutect2 -R $Ref_hg38 -I ${tumor_bam} -I ${normal_bam} -normal NORMAL --germline-resource $germline_resource --panel-of-normals $pon -O Mutect2/${chr}.vcf --native-pair-hmm-threads ${params.threads} -L ${chr}
    """
}

// Process to run HaplotypeCaller for each chromosome
process runHaplotypeCaller {
    input:
    val chr

    output:
    path "HaplotypeCaller/${chr}.vcf"

    script:
    """
    mkdir -p HaplotypeCaller
    source /home/sombuddha/anaconda3/etc/profile.d/conda.sh
    conda activate wgs
    gatk --java-options "-Xmx5g" HaplotypeCaller -R $Ref_hg38 -I ${tumor_bam} --dbsnp $dbsnp -O HaplotypeCaller/${chr}.vcf --native-pair-hmm-threads ${params.threads} -L ${chr}
    """
}

// Workflow definition
workflow {
    chromosomes_ch = Channel.of(*chromosomes)
    mutect2 = runMutect2(chromosomes_ch)
    haplotypecaller = runHaplotypeCaller(chromosomes_ch)
}


