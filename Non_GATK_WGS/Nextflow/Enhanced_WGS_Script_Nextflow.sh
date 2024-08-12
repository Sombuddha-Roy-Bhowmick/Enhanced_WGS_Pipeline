input1="$1"     # Taking input of normal folder name as argument
input2="$2"     # Taking input of tumor folder name as argument

nextflow run Fastq_To_Bam.nf --input1 $1 --input2 $2 --threads 100

./Folder_Arrangement.sh

nextflow run Variant_Calling.nf --threads 100

./Folder_Arrangement_1.sh

rm -rf work/


