#!/bin/bash

# Define base directory as the current working directory
base_dir="$(pwd)"

# Define working directory
work_dir="work"

# Define destination directories
normal_fastp_dir="${base_dir}/Normal_Fastp"
tumor_fastp_dir="${base_dir}/Tumor_Fastp"
tumor_alignment_dir="${base_dir}/Tumor_Alignment"
normal_alignment_dir="${base_dir}/Normal_Alignment"
normal_bqsr_dir="${base_dir}/Normal_BQSR"
tumor_bqsr_dir="${base_dir}/Tumor_BQSR"

# Create the destination directories if they do not exist
mkdir -p "$normal_fastp_dir"
mkdir -p "$tumor_fastp_dir"
mkdir -p "$tumor_alignment_dir"
mkdir -p "$normal_alignment_dir"
mkdir -p "$normal_bqsr_dir"
mkdir -p "$tumor_bqsr_dir"

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

# Move Fastp files for Normal samples
move_files_if_all_present "fastp_normal.html" "fastp_normal.json" "lane1_read1.fastq" "lane1_read2.fastq"
mv "${base_dir}/fastp_normal.html" "${normal_fastp_dir}/"
mv "${base_dir}/fastp_normal.json" "${normal_fastp_dir}/"
mv "${base_dir}/lane1_read1.fastq" "${normal_fastp_dir}/"
mv "${base_dir}/lane1_read2.fastq" "${normal_fastp_dir}/"

# Move Fastp files for Tumor samples
move_files_if_all_present "fastp_tumor.html" "fastp_tumor.json" "lane1_read1.fastq" "lane1_read2.fastq"
mv "${base_dir}/fastp_tumor.html" "${tumor_fastp_dir}/"
mv "${base_dir}/fastp_tumor.json" "${tumor_fastp_dir}/"
mv "${base_dir}/lane1_read1.fastq" "${tumor_fastp_dir}/"
mv "${base_dir}/lane1_read2.fastq" "${tumor_fastp_dir}/"

# Move BAM files from tumor alignment
find "$work_dir" -type d -path "*/tumor_output/Aligned/Projects/default/default" | while read -r dir; do
    echo "Moving BAM files from directory: $dir"
    find "$dir" -type f -name "*.bam*" -exec mv {} "$tumor_alignment_dir/" \;
done

# Move BAM files from normal alignment
find "$work_dir" -type d -path "*/normal_output/Aligned/Projects/default/default" | while read -r dir; do
    echo "Moving BAM files from directory: $dir"
    find "$dir" -type f -name "*.bam*" -exec mv {} "$normal_alignment_dir/" \;
done

# Move BQSR files for Normal samples
move_files_if_all_present "normal_filtered.bam.bai"
mv "${base_dir}/normal_filtered.bam.bai" "${normal_bqsr_dir}/"

find "$work_dir" -type d -path "*/elprep_normal_output" | while read -r dir; do
    echo "Moving BAM files from directory: $dir"
    find "$dir" -type f -name "*" -exec mv {} "$normal_bqsr_dir/" \;
done

# Move BQSR files for Tumor samples
move_files_if_all_present "tumor_filtered.bam.bai"
mv "${base_dir}/tumor_filtered.bam.bai" "${tumor_bqsr_dir}/"

find "$work_dir" -type d -path "*/elprep_tumor_output" | while read -r dir; do
    echo "Moving BAM files from directory: $dir"
    find "$dir" -type f -name "*" -exec mv {} "$tumor_bqsr_dir/" \;
done


echo "Files have been moved to their respective directories."

rm -rf *bam*
rm -rf *.fastq
rm -rf fastp*
rm -rf output.recal

