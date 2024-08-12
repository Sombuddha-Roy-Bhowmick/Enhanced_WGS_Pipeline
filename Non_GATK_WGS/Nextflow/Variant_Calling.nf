#!/usr/bin/env nextflow

// Define parameters with default values
params.threads = params.threads ?: 'default_value'
def ref_hg38 = "Homo_sapiens_assembly38.fasta"
params.workdir = "$PWD"

// Channels for normal and tumor BAM files
def normal_bam = "${params.workdir}/Normal_BQSR/normal_filtered.bam"
def tumor_bam = "${params.workdir}/Tumor_BQSR/tumor_filtered.bam"

// Process to configure Manta workflow
process config_manta {
    input:
    path normal_bam
    path tumor_bam

    output:
    path "MantaWorkflow"

    script:
    """
    mkdir -p MantaWorkflow
    source /home/sombuddha/anaconda3/etc/profile.d/conda.sh
    conda activate py2
    python2.7 /data/softwares/manta-1.6.0.centos6_x86_64/bin/configManta.py --normalBam ${normal_bam} --tumourBam ${tumor_bam} --referenceFasta ${ref_hg38} --runDir MantaWorkflow
    """
}

// Process to run Manta workflow
process run_manta {
    input:
    path manta_config

    output:
    path "MantaWorkflow/results/variants"

    script:
    """
    source /home/sombuddha/anaconda3/etc/profile.d/conda.sh
    conda activate py2
    python2.7  MantaWorkflow/runWorkflow.py -m local -j ${params.threads}
    """
}

// Process to configure Strelka2 Somatic workflow
process config_strelka {
    input:
    path normal_bam
    path tumor_bam
    path manta_out

    output:
    path "StrelkaSomaticWorkflow"

    script:
    """
    source /home/sombuddha/anaconda3/etc/profile.d/conda.sh
    conda activate py2
    python2.7 /data/softwares/strelka-2.9.10.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py --normalBam ${normal_bam} --tumourBam ${tumor_bam} --referenceFasta ${ref_hg38} --indelCandidates ${manta_out}/candidateSmallIndels.vcf.gz  --runDir StrelkaSomaticWorkflow
    """
}

// Process to run Strelka2 Somatic workflow
process run_strelka {
    input:
    path strelka_config

    output:
    path "StrelkaSomaticWorkflow"

    script:
    """
    source /home/sombuddha/anaconda3/etc/profile.d/conda.sh
    conda activate py2
    python2.7  StrelkaSomaticWorkflow/runWorkflow.py -m local -j ${params.threads}
    """
}

// Process to configure Strelka2 Germline workflow
process config_strelka_germline {
    input:
    path normal_bam

    output:
    path "StrelkaGermlineWorkflow"

    script:
    """
    mkdir -p StrelkaGermlineWorkflow
    source /home/tyrone/anaconda3/etc/profile.d/conda.sh
    conda activate py2
    python2.7 /data/softwares/strelka-2.9.10.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py --bam ${normal_bam} --referenceFasta ${ref_hg38} --runDir StrelkaGermlineWorkflow
    """
}

// Process to run Strelka2 Germline workflow
process run_strelka_germline {
    input:
    path strelka_config_germline

    output:
    path "StrelkaGermlineWorkflow"

    script:
    """
    source /home/sombuddha/anaconda3/etc/profile.d/conda.sh
    conda activate py2
    python2.7  StrelkaGermlineWorkflow/runWorkflow.py -m local -j ${params.threads}
    """
}



// Define the workflow
workflow {
    manta_config = config_manta(normal_bam, tumor_bam)
    manta_result = run_manta(manta_config)
    
    strelka_config = config_strelka(normal_bam, tumor_bam, manta_result)
    strelka_result = run_strelka(strelka_config)

    strelka_config_germline = config_strelka_germline(normal_bam)
    strelka_result_germline = run_strelka_germline(strelka_config_germline)
}

