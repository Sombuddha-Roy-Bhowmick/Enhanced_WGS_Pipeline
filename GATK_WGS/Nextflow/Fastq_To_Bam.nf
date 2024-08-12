#!/usr/bin/env nextflow

// Define parameters with default values
params.input1 = params.input1 ?: 'default_value'
params.input2 = params.input2 ?: 'default_value'
params.threads = params.threads ?: 'default_value'
params.workdir = "$PWD"

def raw_normal = "${params.workdir}/${params.input1}"
def raw_tumor = "${params.workdir}/${params.input2}"
def Ref_hg38 = "Homo_sapiens_assembly38.fasta"
def known_sites_1 = "Homo_sapiens_assembly38.dbsnp138.vcf"
def known_sites_2 = "Homo_sapiens_assembly38.known_indels.vcf.gz"



// Create channels for normal and tumor samples
Channel.fromPath("${raw_normal}/*R1_001.fastq.gz").set { normal_r1_files }
Channel.fromPath("${raw_normal}/*R2_001.fastq.gz").set { normal_r2_files }
Channel.fromPath("${raw_tumor}/*R1_001.fastq.gz").set { tumor_r1_files }
Channel.fromPath("${raw_tumor}/*R2_001.fastq.gz").set { tumor_r2_files }

// Process to run fastp on normal samples
process fastp_normal {
    
    input:
    path r1_files
    path r2_files
    
    output:
    tuple path("normal_trimmed_R1.fastq"), path("normal_trimmed_R2.fastq")
    
    script:
    """
    fastp --detect_adapter_for_pe --overrepresentation_analysis --correction --cut_right --thread ${params.threads} \
          --html fastp_normal.html \
          --json fastp_normal.json \
          -i ${r1_files} \
          -o normal_trimmed_R1.fastq \
          -I ${r2_files} \
          -O normal_trimmed_R2.fastq
    """
}

process fastp_tumor {
    
    input:
    path r1_files
    path r2_files
    
    output:
    tuple path("tumor_trimmed_R1.fastq"), path("tumor_trimmed_R2.fastq")
    
    script:
    """
    fastp --detect_adapter_for_pe --overrepresentation_analysis --correction --cut_right --thread ${params.threads} \
          --html fastp_tumor.html \
          --json fastp_tumor.json \
          -i ${r1_files} \
          -o tumor_trimmed_R1.fastq \
          -I ${r2_files} \
          -O tumor_trimmed_R2.fastq
    """
}

// Process to run bwa mem on normal samples
process bwa_mem_normal {
    input:
    tuple path(trimmed_r1), path(trimmed_r2)
    
    output:
    path "Normal_Alignment/normal_sorted.bam"
    
    script:
    """
    mkdir -p Normal_Alignment
    source /home/tyrone/anaconda3/etc/profile.d/conda.sh
    conda activate wgs
    bwa mem -t ${params.threads} -M -R "@RG\\tID:NORMAL\\tPL:ILLUMINA\\tLB:TruSeq\\tSM:NORMAL\\tPI:200" $Ref_hg38 ${trimmed_r1} ${trimmed_r2} | samtools sort -@ ${params.threads} -o Normal_Alignment/normal_sorted.bam
    """
}

// Process to run bwa mem on tumor samples
process bwa_mem_tumor {
    input:
    tuple path(trimmed_r1), path(trimmed_r2)
    
    output:
    path "Tumor_Alignment/tumor_sorted.bam"
    
    script:
    """
    mkdir -p Tumor_Alignment
    source /home/tyrone/anaconda3/etc/profile.d/conda.sh
    conda activate wgs
    bwa mem -t ${params.threads} -M -R "@RG\\tID:TUMOR\\tPL:ILLUMINA\\tLB:TruSeq\\tSM:TUMOR\\tPI:200" $Ref_hg38 ${trimmed_r1} ${trimmed_r2} | samtools sort -@ ${params.threads} -o Tumor_Alignment/tumor_sorted.bam
    """
}

// Process to run MarkDuplicates on normal samples
process markduplicates_normal {
    input:
    path sorted_bam
    
    output:
    path "Normal_Markduplicate/Normal_dedup.bam"
    
    script:
    """
    mkdir -p Normal_Markduplicate
    source /home/tyrone/anaconda3/etc/profile.d/conda.sh
    conda activate wgs
    gatk --java-options "-Xms100g" MarkDuplicates -I ${sorted_bam} -O Normal_Markduplicate/Normal_dedup.bam -M Normal_Markduplicate/normal_metrics.txt --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 -R $Ref_hg38
    """
}

// Process to run MarkDuplicates on tumor samples
process markduplicates_tumor {
    input:
    path sorted_bam
    
    output:
    path "Tumor_Markduplicate/Tumor_dedup.bam"
    
    script:
    """
    mkdir -p Tumor_Markduplicate
    source /home/tyrone/anaconda3/etc/profile.d/conda.sh
    conda activate wgs
    gatk --java-options "-Xms100g" MarkDuplicates -I ${sorted_bam} -O Tumor_Markduplicate/Tumor_dedup.bam -M Tumor_Markduplicate/tumor_metrics.txt --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 -R $Ref_hg38
    """
}

// Process to run BaseRecalibrator on normal samples
process baserecalibrator_normal {
    input:
    path dedup_bam
    
    output:
    path "Normal_BQSR/normal_recal.table"
    
    script:
    """
    mkdir -p Normal_BQSR
    source /home/tyrone/anaconda3/etc/profile.d/conda.sh
    conda activate wgs
    gatk --java-options "-Xmx100g" BaseRecalibrator --input ${dedup_bam} --reference $Ref_hg38 --known-sites $known_sites_1 --known-sites $known_sites_2 --output Normal_BQSR/normal_recal.table
    """
}

// Process to run BaseRecalibrator on tumor samples
process baserecalibrator_tumor {
    input:
    path dedup_bam
    
    output:
    path "Tumor_BQSR/tumor_recal.table"
    
    script:
    """
    mkdir -p Tumor_BQSR
    source /home/tyrone/anaconda3/etc/profile.d/conda.sh
    conda activate wgs
    gatk --java-options "-Xmx100g" BaseRecalibrator --input ${dedup_bam} --reference $Ref_hg38 --known-sites $known_sites_1 --known-sites $known_sites_2 --output Tumor_BQSR/tumor_recal.table
    """
}

// Process to run ApplyBQSR on normal samples
process applybqsr_normal {
    input:
    tuple path(dedup_bam), path(recal_table)
    
    output:
    path "BQSR_Normal/normal_bqsr.bam"
    
    script:
    """
    mkdir -p BQSR_Normal
    source /home/tyrone/anaconda3/etc/profile.d/conda.sh
    conda activate wgs
    gatk --java-options "-Xmx100g" ApplyBQSR -R $Ref_hg38 -I ${dedup_bam} -bqsr ${recal_table} -O BQSR_Normal/normal_bqsr.bam
    """
}

// Process to run ApplyBQSR on tumor samples
process applybqsr_tumor {
    input:
    tuple path(dedup_bam), path(recal_table)
    
    output:
    path "BQSR_Tumor/tumor_bqsr.bam"
    
    script:
    """
    mkdir -p BQSR_Tumor
    source /home/tyrone/anaconda3/etc/profile.d/conda.sh
    conda activate wgs
    gatk --java-options "-Xmx100g" ApplyBQSR -R $Ref_hg38 -I ${dedup_bam} -bqsr ${recal_table} -O BQSR_Tumor/tumor_bqsr.bam
    """
}

// Define workflow
 
workflow {
    
    normal_r1 = Channel.fromPath("${raw_normal}/*R1_001.fastq.gz")
    normal_r2 = Channel.fromPath("${raw_normal}/*R2_001.fastq.gz")
    tumor_r1 = Channel.fromPath("${raw_tumor}/*R1_001.fastq.gz")
    tumor_r2 = Channel.fromPath("${raw_tumor}/*R2_001.fastq.gz")

    fastp_normal_results = fastp_normal(normal_r1, normal_r2)
    fastp_tumor_results = fastp_tumor(tumor_r1, tumor_r2)

    bwa_mem_normal_results = bwa_mem_normal(fastp_normal_results)
    bwa_mem_tumor_results = bwa_mem_tumor(fastp_tumor_results)

    markduplicates_normal_results = markduplicates_normal(bwa_mem_normal_results)
    markduplicates_tumor_results = markduplicates_tumor(bwa_mem_tumor_results)
   
    baserecalibrator_normal_results = baserecalibrator_normal(markduplicates_normal_results)
    baserecalibrator_tumor_results = baserecalibrator_tumor(markduplicates_tumor_results)
  
    // Combine markduplicates_normal_results and baserecalibrator_normal_results into a tuple
    normal_bqsr_inputs = markduplicates_normal_results.combine(baserecalibrator_normal_results)
    tumor_bqsr_inputs = markduplicates_tumor_results.combine(baserecalibrator_tumor_results)

    applybqsr_normal_results = applybqsr_normal(normal_bqsr_inputs)
    applybqsr_tumor_results = applybqsr_tumor(tumor_bqsr_inputs)   
}

