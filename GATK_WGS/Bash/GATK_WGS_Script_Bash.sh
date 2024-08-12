input1="$1"     # Taking input of normal folder name as argument
input2="$2"     # Taking input of tumor folder name as argument
workdir=$(pwd)  # storing path of Working directory
raw_normal="$workdir/$input1"   # Assigning path of normal sample
raw_tumor="$workdir/$input2"    # Assigning path of tumor sample 
out="$workdir/Output_hg38_${input1}_Vs_${input2}"       # Assigning path of output folder
Ref_hg38="Homo_sapiens_assembly38.fasta" # Reference Genome for hg38
known_sites_1="Homo_sapiens_assembly38.dbsnp138.vcf" # Known Sites - 1
known_sites_2="Homo_sapiens_assembly38.known_indels.vcf.gz" # Known Sites - 2
threads="100"
half_threads="50"
germline_resource="af-only-gnomad.hg38.vcf.gz" # Germline Resource 
pon="1000g_pon.hg38.vcf.gz"

cd $input1
ls *R1* > ../list_${input1}

cd $workdir
input3="list_${input1}"   # Creating list file for normal sample

cd $input2
ls *R1* > ../list_${input2}

cd $workdir
input4="list_${input2}"  # Creating list file for tumor sample


mkdir -p $out/$input1/fastp_output
mkdir -p $out/$input1/Alignment
mkdir -p $out/$input1/BQSR
mkdir -p $out/$input2/fastp_output
mkdir -p $out/$input2/Alignment
mkdir -p $out/$input2/BQSR
mkdir -p $out/variantcaller/Mutect2
mkdir -p $out/variantcaller/HaplotypeCaller






eval "$(conda shell.bash hook)"         # Setting bash for conda environment
conda activate wgs             # Activating conda environment

#Fastp for Normal


echo "fastp --detect_adapter_for_pe --overrepresentation_analysis --correction --cut_right --thread $threads --html $out/$input1/fastp_output/${input1}.fastp.html --json $out/$input1/fastp_output/${input1}.fastp.json -i $raw_normal/${input1}*R1_001.fastq.gz  -o $out/$input1/fastp_output/${input1}_trimmed_R1.fastq -I $raw_normal/${input1}*R2_001.fastq.gz  -O $out/$input1/fastp_output/${input1}_trimmed_R2.fastq" >> $out/parallel_fastp


#Fastp for Tumor

echo "fastp --detect_adapter_for_pe --overrepresentation_analysis --correction --cut_right --thread $threads  --html $out/$input2/fastp_output/${input2}.fastp.html --json $out/$input2/fastp_output/${input2}.fastp.json -i $raw_tumor/${input2}*R1_001.fastq.gz  -o $out/$input2/fastp_output/${input2}_trimmed_R1.fastq -I $raw_tumor/${input2}*R2_001.fastq.gz  -O $out/$input2/fastp_output/${input2}_trimmed_R2.fastq" >> $out/parallel_fastp



parallel -j 2 < $out/parallel_fastp # Running fastp  in parallel

#Mapping for Normal 

bwa mem -t $threads  -M -R "@RG\tID:NORMAL\tPL:ILLUMINA\tLB:TruSeq\tSM:NORMAL\tPI:200" $Ref_hg38 $out/$input1/fastp_output/${input1}_trimmed_R1.fastq  $out/$input1/fastp_output/${input1}_trimmed_R2.fastq | samtools sort -@ $threads  -o $out/$input1/Alignment/${input1}_sorted.bam

#Mapping for Tumor

bwa mem -t $threads  -M -R "@RG\tID:TUMOR\tPL:ILLUMINA\tLB:TruSeq\tSM:TUMOR\tPI:200" $Ref_hg38 $out/$input2/fastp_output/${input2}_trimmed_R1.fastq  $out/$input2/fastp_output/${input2}_trimmed_R2.fastq | samtools sort -@ $threads  -o $out/$input2/Alignment/${input2}_sorted.bam



#Markduplicate for Normal

echo "gatk  --java-options "-Xms100g"  MarkDuplicates  -I $out/$input1/Alignment/${input1}_sorted.bam -O $out/$input1/Alignment/${input1}_dedup.bam -M $out/$input1/Alignment/${input1}_metrics.txt --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 -R $Ref_hg38" >>  $out/parallel_markdup

#Markduplicate for Tumor

echo "gatk  --java-options "-Xms100g"  MarkDuplicates -I $out/$input2/Alignment/${input2}_sorted.bam -O $out/$input2/Alignment/${input2}_dedup.bam -M $out/$input2/Alignment/${input2}_metrics.txt --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 -R $Ref_hg38" >>  $out/parallel_markdup

parallel -j 2 < $out/parallel_markdup  # Running markduplicates  in parallel


#Base Recalibration for Normal

echo "gatk  --java-options "-Xmx100g" BaseRecalibrator --input $out/$input1/Alignment/${input1}_dedup.bam --reference $Ref_hg38 --known-sites $known_sites_1 --known-sites $known_sites_2 --output $out/$input1/BQSR/${input1}_recal.table" >> $out/parallel_baserecal

#Base Recalibration for Tumor

echo "gatk  --java-options "-Xmx100g" BaseRecalibrator --input $out/$input2/Alignment/${input2}_dedup.bam --reference $Ref_hg38 --known-sites $known_sites_1 --known-sites $known_sites_2 --output $out/$input2/BQSR/${input2}_recal.table" >> $out/parallel_baserecal

parallel -j 2 < $out/parallel_baserecal # Running baserecal  in parallel

#Apply BQSR for Normal

echo "gatk --java-options "-Xmx100g" ApplyBQSR -R $Ref_hg38 -I $out/$input1/Alignment/${input1}_dedup.bam  -bqsr $out/$input1/BQSR/${input1}_recal.table -O $out/$input1/BQSR/${input1}_bqsr.bam" >> $out/parallel_applybqsr

#Apply BQSR for Tumor

echo "gatk --java-options "-Xmx100g" ApplyBQSR -R $Ref_hg38 -I $out/$input2/Alignment/${input2}_dedup.bam  -bqsr $out/$input2/BQSR/${input2}_recal.table -O $out/$input2/BQSR/${input2}_bqsr.bam" >> $out/parallel_applybqsr

parallel -j 2 < $out/parallel_applybqsr # Running applybqsr  in parallel


#Somatic Vaiant Calling

pids_mutect=
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY
	do
		gatk --java-options "-Xmx5g" Mutect2 -R $Ref_hg38  -I $out/$input2/BQSR/${input2}_bqsr.bam -I $out/$input1/BQSR/${input1}_bqsr.bam  -normal NORMAL  --germline-resource $germline_resource --panel-of-normals $pon -O $out/variantcaller/Mutect2/$chr.vcf --native-pair-hmm-threads 10 -L $chr &
	pids_mutect+=" $!"
done;

wait $pids_mutect
wait


#Germline Variant Calling

pids=
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY
      do
      gatk --java-options "-Xmx5g" HaplotypeCaller -R $Ref_hg38  -I $out/$input2/BQSR/${input2}_bqsr.bam --dbsnp $known_sites_1  -O $out/variantcaller/HaplotypeCaller/$chr.vcf --native-pair-hmm-threads 10 -L $chr &
      pids+=" $!"
    done;

wait $pids
wait


#Somatic VCF Merge
gatk MergeVcfs I=$out/variantcaller/Mutect2/chr1.vcf I=$out/variantcaller/Mutect2/chr2.vcf I=$out/variantcaller/Mutect2/chr3.vcf I=$out/variantcaller/Mutect2/chr4.vcf I=$out/variantcaller/Mutect2/chr5.vcf I=$out/variantcaller/Mutect2/chr6.vcf I=$out/variantcaller/Mutect2/chr7.vcf I=$out/variantcaller/Mutect2/chr8.vcf I=$out/variantcaller/Mutect2/chr9.vcf I=$out/variantcaller/Mutect2/chr10.vcf I=$out/variantcaller/Mutect2/chr11.vcf I=$out/variantcaller/Mutect2/chr12.vcf I=$out/variantcaller/Mutect2/chr13.vcf I=$out/variantcaller/Mutect2/chr14.vcf I=$out/variantcaller/Mutect2/chr15.vcf I=$out/variantcaller/Mutect2/chr16.vcf I=$out/variantcaller/Mutect2/chr17.vcf I=$out/variantcaller/Mutect2/chr18.vcf I=$out/variantcaller/Mutect2/chr19.vcf I=$out/variantcaller/Mutect2/chr20.vcf I=$out/variantcaller/Mutect2/chr21.vcf I=$out/variantcaller/Mutect2/chr22.vcf I=$out/variantcaller/Mutect2/chrX.vcf I=$out/variantcaller/Mutect2/chrY.vcf O=$out/variantcaller/Mutect2/Mutect2.vcf

#Germline VCF Merge

gatk MergeVcfs I=$out/variantcaller/HaplotypeCaller/chr1.vcf I=$out/variantcaller/HaplotypeCaller/chr2.vcf I=$out/variantcaller/HaplotypeCaller/chr3.vcf I=$out/variantcaller/HaplotypeCaller/chr4.vcf I=$out/variantcaller/HaplotypeCaller/chr5.vcf I=$out/variantcaller/HaplotypeCaller/chr6.vcf I=$out/variantcaller/HaplotypeCaller/chr7.vcf I=$out/variantcaller/HaplotypeCaller/chr8.vcf I=$out/variantcaller/HaplotypeCaller/chr9.vcf I=$out/variantcaller/HaplotypeCaller/chr10.vcf I=$out/variantcaller/HaplotypeCaller/chr11.vcf I=$out/variantcaller/HaplotypeCaller/chr12.vcf I=$out/variantcaller/HaplotypeCaller/chr13.vcf I=$out/variantcaller/HaplotypeCaller/chr14.vcf I=$out/variantcaller/HaplotypeCaller/chr15.vcf I=$out/variantcaller/HaplotypeCaller/chr16.vcf I=$out/variantcaller/HaplotypeCaller/chr17.vcf I=$out/variantcaller/HaplotypeCaller/chr18.vcf I=$out/variantcaller/HaplotypeCaller/chr19.vcf I=$out/variantcaller/HaplotypeCaller/chr20.vcf I=$out/variantcaller/HaplotypeCaller/chr21.vcf I=$out/variantcaller/HaplotypeCaller/chr22.vcf I=$out/variantcaller/HaplotypeCaller/chrX.vcf I=$out/variantcaller/HaplotypeCaller/chrY.vcf O=$out/variantcaller/HaplotypeCaller/HaplotypeCaller.vcf












