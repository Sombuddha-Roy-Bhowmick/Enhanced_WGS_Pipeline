# Enhanced WGS Pipeline

**The Enhanced Whole Genome Sequencing Pipeline is ~4x faster than the standard GATK WGS pipeline.**

The Enhanced Non-GATK pipeline uses various open source tools:
1. **fastp** (https://github.com/OpenGene/fastp) - QC, Adapter Removal & Trimming
2. **Isaac aligner** (https://github.com/Illumina/Isaac4) - Mapping/Alignment & Marking Duplicate Reads
3. **elPrep** (https://github.com/ExaScience/elprep) - Base Quality Score Recalibration & InDel Realignment
4. **Strelka2** (https://github.com/Illumina/strelka) - Somatic & Germline Variant Caller

The GATK pipeline uses various open source tools:
1. **fastp** (https://github.com/OpenGene/fastp) - QC, Adapter Removal & Trimming
2. **GATK** (https://github.com/broadinstitute/gatk) - Mapping/Alignment, Marking Duplicate Reads, Base Quality Score Recalibration and InDel Realignment, Somatic (Mutect2) & Germline (HaplotypeCaller) Variant Caller

Both the pipelines are avilable in this repository as **Shell** scripts as well as **Nextflow** scripts.

**For a 80x-70x Tumor-MatchNormal sample, the Enhanced pipleine completed the analysis in 6 hours 26 minutes, whereas the GATK pipeline took 25 hours 32 minutes for the same pair of samples.**

The output of the two pipelines were found to be concordant. The pipelines were tested on the cell line data - HCC1395 (available from NCBI-SRA) and both the pipelines could capture the same number of true somatic mutations (information obtained from COSMIC) at different depths.

The Enhanced pipeline turned out to be much faster than the conventional GATK pipeline.




