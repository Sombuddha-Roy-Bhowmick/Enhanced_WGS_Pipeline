#!/usr/bin/env nextflow

// Define parameters with default values
params.input1 = params.input1 ?: 'default_value'
params.input2 = params.input2 ?: 'default_value'
params.threads = params.threads ?: 'default_value'
params.workdir = "$PWD"

def raw_normal = "${params.workdir}/${params.input1}"
def raw_tumor = "${params.workdir}/${params.input2}"
def ref_hg38 = "/data/002_genomes/WGS_Reference_hg38/GCA_000001405.29_GRCh38.p14_genomic.fa"
def isaac_hg38 = "/data/002_genomes/WGS_Reference_hg38/IsaacIndex.20240724/sorted-reference.xml"
def elprep_ref = "/data/002_genomes/WGS_Reference_hg38/GCA_000001405.29_GRCh38.p14_genomic.elfasta"
def elprep_known_sites = "/data/002_genomes/WGS_Reference_hg38/Homo_sapiens_assembly_38.known_indels.elsites,/data/002_genomes/WGS_Reference_hg38/Homo_sapiens_assembly38.dbsnp138.elsites"

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
    tuple path("lane1_read1.fastq"), path("lane1_read2.fastq") 

    script:
    """
    fastp --detect_adapter_for_pe --overrepresentation_analysis --correction --cut_right --thread ${params.threads} \
          --html fastp_normal.html \
          --json fastp_normal.json \
          -i ${r1_files} \
          -o lane1_read1.fastq \
          -I ${r2_files} \
          -O lane1_read2.fastq
    """
}

// Process to run fastp on tumor samples
process fastp_tumor {
    input:
    path r1_files
    path r2_files

    output:
    tuple path("lane1_read1.fastq"), path("lane1_read2.fastq")

    script:
    """
    fastp --detect_adapter_for_pe --overrepresentation_analysis --correction --cut_right --thread ${params.threads} \
          --html fastp_tumor.html \
          --json fastp_tumor.json \
          -i ${r1_files} \
          -o lane1_read1.fastq \
          -I ${r2_files} \
          -O lane1_read2.fastq
    """
}

// Process to run isaac-align on normal samples
process isaac_align_normal {
    input:
    tuple path(read1), path(read2)

    output:
    path "normal_output/Aligned/Projects/default/default/"

    script:
    """
    mkdir -p normal_output/Aligned/Projects/default/default
    mv ${read1} ${read2} normal_output/
    source /home/tyrone/anaconda3/etc/profile.d/conda.sh
    conda activate wgs
    /data/007_Softwares/Isaac4/Isaac-build/bin/isaac-align -r ${isaac_hg38} -b normal_output -f fastq --use-bases-mask y150,y150 -m200 -j ${params.threads} --keep-duplicates 1 --mark-duplicates 1 -o normal_output/Aligned/
    """
}

// Process to run isaac-align on tumor samples
process isaac_align_tumor {
    input:
    tuple path(read1), path(read2)

    output:
    path "tumor_output/Aligned/Projects/default/default/"

    script:
    """
    mkdir -p tumor_output/Aligned/Projects/default/default
    mv ${read1} ${read2} tumor_output/
    source /home/tyrone/anaconda3/etc/profile.d/conda.sh
    conda activate wgs
    /data/007_Softwares/Isaac4/Isaac-build/bin/isaac-align -r ${isaac_hg38} -b tumor_output -f fastq --use-bases-mask y150,y150 -m200 -j ${params.threads} --keep-duplicates 1 --mark-duplicates 1 -o tumor_output/Aligned/
    ls -l tumor_output/Aligned/Projects/default/default
    """
}

// Process to run elprep filter on normal samples
process elprep_normal {
    input:
    path bam_input

    output:
    tuple path("elprep_normal_output/normal_filtered.bam"), path("elprep_normal_output/output.recal")

    script:
    """
    mkdir -p elprep_normal_output
    source /home/tyrone/anaconda3/etc/profile.d/conda.sh
    conda activate wgs
    elprep sfm  ${bam_input}/*bam  elprep_normal_output/normal_filtered.bam  --bqsr elprep_normal_output/output.recal --reference ${elprep_ref} --known-sites ${elprep_known_sites} --nr-of-threads ${params.threads}
    """
}

// Process to run elprep filter on tumor samples
process elprep_tumor {
    input:
    path bam_input

    output:
    tuple path("elprep_tumor_output/tumor_filtered.bam"), path("elprep_tumor_output/output.recal")

    script:
    """
    mkdir -p elprep_tumor_output
    source /home/tyrone/anaconda3/etc/profile.d/conda.sh
    conda activate wgs
    elprep sfm  ${bam_input}/*bam  elprep_tumor_output/tumor_filtered.bam  --bqsr elprep_tumor_output/output.recal --reference ${elprep_ref} --known-sites ${elprep_known_sites} --nr-of-threads ${params.threads}
    """
}

// Process to index the BAM files for normal samples using samtools
process index_normal {
    input:
    tuple path(bam_file), path(_)

    output:
    path "${bam_file}.bai"

    script:
    """
    mkdir -p elprep_normal_output
    source /home/tyrone/anaconda3/etc/profile.d/conda.sh
    conda activate wgs
    samtools index -@ ${params.threads} ${bam_file.toString()}
    """
}

// Process to index the BAM files for tumor samples using samtools
process index_tumor {
    input:
     tuple path(bam_file), path(_)

    output:
    path "${bam_file}.bai"

    script:
    """
    mkdir -p elprep_tumor_output
    source /home/tyrone/anaconda3/etc/profile.d/conda.sh
    conda activate wgs
    samtools index  -@ ${params.threads} ${bam_file.toString()}
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

    isaac_normal_bam = isaac_align_normal(fastp_normal_results)
    isaac_tumor_bam = isaac_align_tumor(fastp_tumor_results)

    // Run elprep for normal and tumor samples
    elprep_normal_output = elprep_normal(isaac_normal_bam)
    elprep_tumor_output = elprep_tumor(isaac_tumor_bam)
    
    index_normal(elprep_normal_output)
    index_tumor(elprep_tumor_output)
}
