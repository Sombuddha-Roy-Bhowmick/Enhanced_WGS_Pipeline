#!/bin/bash

# Define base directory as the current working directory
base_dir="$(pwd)"

# Define working directory
work_dir="work"

# Define destination directories
Mutect2_dir="${base_dir}/Mutect2"
HaplotypeCaller_dir="${base_dir}/HaplotypeCaller"

# Create the destination directories if they do not exist
mkdir -p "$Mutect2_dir"
mkdir -p "$HaplotypeCaller_dir"

# Function to move files if all required files are present
move_files_if_all_present() {
    local target_files=("$@")
    local target_dir

    for target_dir in $(find "$work_dir" -type d); do
        if [[ -d "$target_dir" ]]; then
            local all_present=true
            for file in "${target_files[@]}"; do
                if [[ ! -f "$target_dir/$file" ]]; then
                    all_present=false
                    break
                fi
            done

            if $all_present; then
                echo "Moving files from directory: $target_dir"
                for file in "${target_files[@]}"; do
                    mv "$target_dir/$file" "$base_dir"
                done
            fi
        fi
    done
}

# Function to resolve symbolic links and move files
move_files() {
    local source_dir="$1"
    local dest_dir="$2"
    local files=("$@")
    files=("${files[@]:2}")

    mkdir -p "$dest_dir"

    for file in "${files[@]}"; do
        # Resolve symbolic link or use file directly
        if [ -L "$source_dir/$file" ]; then
            actual_file=$(readlink -f "$source_dir/$file")
        else
            actual_file="$source_dir/$file"
        fi

        if [ -f "$actual_file" ]; then
            echo "Moving file: $actual_file to $dest_dir"
            mv "$actual_file" "$dest_dir/"
        fi
    done
}


# Move VCF files from Mutect2
find "$work_dir" -type d -path "*/Mutect2" | while read -r dir; do
    echo "Moving VCF files from directory: $dir"
    find "$dir" -type f -name "*vcf*" -exec mv {} "$Mutect2_dir/" \;
done

# Move VCF files from HaplotypeCaller
find "$work_dir" -type d -path "*/HaplotypeCaller" | while read -r dir; do
    echo "Moving VCF files from directory: $dir"
    find "$dir" -type f -name "*vcf*" -exec mv {} "$HaplotypeCaller_dir/" \;
done


echo "Files have been moved to their respective directories."

#rm -rf *bam*
#rm -rf *.fastq
#rm -rf fastp*
#rm -rf output.recal


eval "$(conda shell.bash hook)"         # Setting bash for conda environment
conda activate wgs             # Activating conda environment


#Merge Mutect2 VCF

gatk MergeVcfs I=Mutect2/chr1.vcf I=Mutect2/chr2.vcf I=Mutect2/chr3.vcf I=Mutect2/chr4.vcf I=Mutect2/chr5.vcf I=Mutect2/chr6.vcf I=Mutect2/chr7.vcf I=Mutect2/chr8.vcf I=Mutect2/chr9.vcf I=Mutect2/chr10.vcf I=Mutect2/chr11.vcf I=Mutect2/chr12.vcf I=Mutect2/chr13.vcf I=Mutect2/chr14.vcf I=Mutect2/chr15.vcf I=Mutect2/chr16.vcf I=Mutect2/chr17.vcf I=Mutect2/chr18.vcf I=Mutect2/chr19.vcf I=Mutect2/chr20.vcf I=Mutect2/chr21.vcf I=Mutect2/chr22.vcf I=Mutect2/chrX.vcf I=Mutect2/chrY.vcf O=Mutect2/Mutect2.vcf

#Merge HaplotypeCaller VCF
gatk MergeVcfs I=HaplotypeCaller/chr1.vcf I=HaplotypeCaller/chr2.vcf I=HaplotypeCaller/chr3.vcf I=HaplotypeCaller/chr4.vcf I=HaplotypeCaller/chr5.vcf I=HaplotypeCaller/chr6.vcf I=HaplotypeCaller/chr7.vcf I=HaplotypeCaller/chr8.vcf I=HaplotypeCaller/chr9.vcf I=HaplotypeCaller/chr10.vcf I=HaplotypeCaller/chr11.vcf I=HaplotypeCaller/chr12.vcf I=HaplotypeCaller/chr13.vcf I=HaplotypeCaller/chr14.vcf I=HaplotypeCaller/chr15.vcf I=HaplotypeCaller/chr16.vcf I=HaplotypeCaller/chr17.vcf I=HaplotypeCaller/chr18.vcf I=HaplotypeCaller/chr19.vcf I=HaplotypeCaller/chr20.vcf I=HaplotypeCaller/chr21.vcf I=HaplotypeCaller/chr22.vcf I=HaplotypeCaller/chrX.vcf I=HaplotypeCaller/chrY.vcf O=HaplotypeCaller/HaplotypeCaller.vcf


