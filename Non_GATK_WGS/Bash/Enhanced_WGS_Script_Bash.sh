input1="$1"	# Taking input of normal folder name as argument
input2="$2"	# Taking input of tumor folder name as argument
workdir=$(pwd)	# storing path of Working directory
raw_normal="$workdir/$input1"	# Assigning path of normal sample
raw_tumor="$workdir/$input2"	# Assigning path of tumor sample 
out="$workdir/Output_hg38_${input1}_Vs_${input2}"	# Assigning path of output folder
Ref_hg38="Homo_sapiens_assembly38.fasta" # Reference Genome for hg38
Isaac_hg38="IsaacIndex.20240801/sorted-reference.xml" # Isaac Reference Genome for hg38
elprep_ref="Homo_sapiens_assembly38.elfasta" # Elprep Reference
elprep_known_sites="Homo_sapiens_assembly_38.known_indels.elsites,Homo_sapiens_assembly38.dbsnp138.elsites" # Known Sites 
threads="100"
half_threads="50"


rm $out/parallel_fastp
rm $out/parallel_bamindex

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
mkdir -p $out/variantCaller/Manta
mkdir -p $out/variantCaller/Strelka2_Somatic
mkdir -p $out/variantCaller/Strelka2_Germline


eval "$(conda shell.bash hook)"         # Setting bash for conda environment
conda activate wgs             # Activating conda environment

while IFS= read -r line			
	do 				
		echo $line
		sm=$(echo $raw_normal/$line | xargs -n 1 basename | cut -f 1 -d"_"); 	# Extracting file name from file 
		echo $sm




	echo "fastp --detect_adapter_for_pe --overrepresentation_analysis --correction --cut_right --thread $threads --html $out/$input1/fastp_output/${sm}.fastp.html --json $out/$input1/fastp_output/${sm}.fastp.json -i $raw_normal/${sm}*_R1_001.fastq.gz  -o $out/$input1/fastp_output/lane1_read1.fastq -I $raw_normal/${sm}*_R2_001.fastq.gz  -O $out/$input1/fastp_output/lane1_read2.fastq" >> $out/parallel_fastp

	done < "$input3"

while IFS= read -r line                 
        do
                echo $line
                sm=$(echo $raw_tumor/$line | xargs -n 1 basename | cut -f 1 -d"_");    # Extracting file name from file 
                echo $sm

	echo "fastp --detect_adapter_for_pe --overrepresentation_analysis --correction --cut_right --thread $threads  --html $out/$input2/fastp_output/${sm}.fastp.html --json $out/$input2/fastp_output/${sm}.fastp.json -i $raw_tumor/${sm}*_R1_001.fastq.gz  -o $out/$input2/fastp_output/lane1_read1.fastq -I $raw_tumor/${sm}*_R2_001.fastq.gz  -O $out/$input2/fastp_output/lane1_read2.fastq" >> $out/parallel_fastp

	done < "$input4"

parallel -j 2 < $out/parallel_fastp # Running fastp  in parallel

#Alignment & Markduplication of Normal Sample

/data/softwares/Isaac4/Isaac-build/bin/isaac-align -r $Isaac_hg38  -b $out/$input1/fastp_output -f fastq --use-bases-mask y150,y150 -m200 -j $threads --keep-duplicates 1 --mark-duplicates 1  -o $out/$input1/Alignment -t $out/$input1/Alignment/Temp

rm -rf $out/$input1/Alignment/Temp

#Alignment & Markduplication of Tumor Sample

/data/softwares/Isaac4/Isaac-build/bin/isaac-align -r $Isaac_hg38 -b $out/$input2/fastp_output -f fastq --use-bases-mask y150,y150 -m200 -j $threads  --keep-duplicates 1 --mark-duplicates 1  -o $out/$input2/Alignment -t $out/$input2/Alignment/Temp

rm -rf $out/$input2/Alignment/Temp

#BQSR - Normal Sample

elprep sfm  $out/$input1/Alignment/Projects/default/default/sorted.bam  $out/$input1/BQSR/${input1}.bqsr.bam  --bqsr $out/$input1/BQSR/${input1}.output.recal --reference $elprep_ref  --known-sites $elprep_known_sites  --nr-of-threads $threads

#BQSR - Tumor Sample

elprep sfm  $out/$input2/Alignment/Projects/default/default/sorted.bam  $out/$input2/BQSR/${input2}.bqsr.bam  --bqsr $out/$input2/BQSR/${input2}.output.recal --reference $elprep_ref  --known-sites $elprep_known_sites  --nr-of-threads $threads

# Bam Indexing - Normal & Tumor Sample

echo "samtools index -@ $half_threads  $out/$input1/BQSR/${input1}.bqsr.bam" >> $out/parallel_bamindex
echo "samtools index -@ $half_threads  $out/$input2/BQSR/${input2}.bqsr.bam" >> $out/parallel_bamindex 

parallel -j 2 < $out/parallel_bamindex


eval "$(conda shell.bash hook)"         # Setting bash for conda environment
conda deactivate             # Deactivating conda environment

eval "$(conda shell.bash hook)"		# Setting bash for conda environment
conda activate py2		# Activating conda environment

#SV Caller

python /data/softwares/manta-1.6.0.centos6_x86_64/bin/configManta.py --normalBam $out/$input1/BQSR/${input1}.bqsr.bam --tumourBam $out/$input2/BQSR/${input2}.bqsr.bam --referenceFasta $Ref_hg38  --runDir $out/variantCaller/Manta

python $out/variantCaller/Manta/runWorkflow.py  -j $threads

#Somatic Variant Caller

python /data/softwares/strelka-2.9.10.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py --normalBam $out/$input1/BQSR/${input1}.bqsr.bam --tumourBam $out/$input2/BQSR/${input2}.bqsr.bam --referenceFasta $Ref_hg38 --indelCandidates $out/variantCaller/Manta/results/variants/candidateSmallIndels.vcf.gz  --runDir $out/variantCaller/Strelka2_Somatic 

python $out/variantCaller/Strelka2_Somatic/runWorkflow.py -m local -j $threads

#Germline Variant Caller

python /data/softwares/strelka-2.9.10.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py --bam $out/$input2/BQSR/${input1}.bqsr.bam --referenceFasta $Ref_hg38  --runDir $out/variantCaller/Strelka2_Germline

python $out/variantCaller/Strelka2_Germline/runWorkflow.py -m local -j $threads


eval "$(conda shell.bash hook)"         # Setting bash for conda environment
conda deactivate             # Deactivating conda environment

