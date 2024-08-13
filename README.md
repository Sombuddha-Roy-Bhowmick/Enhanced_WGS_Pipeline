# Enhanced_WGS_Pipeline

The Enhanced Whole Genome Sequencing Pipeline is ~4x faster than the standard GATK WGS pipeline. 

The Enhanced Non-GATK pipeline uses various open source tools:
1. fastp (https://github.com/OpenGene/fastp) - QC, Adapter Removal & Trimming
2. Isaac aligner (https://github.com/Illumina/Isaac4) - Mapping/Alignment & Marking Duplicate Reads
3. elPrep (https://github.com/ExaScience/elprep) - Base Quality Score Recalibration & InDel Realignment
4. Strelka2 (https://github.com/Illumina/strelka) - Somatic & Germline Variant Caller

The GATK pipeline uses:
1. fastp (https://github.com/OpenGene/fastp) - QC, Adapter Removal & Trimming
2. GATK (https://github.com/broadinstitute/gatk)


