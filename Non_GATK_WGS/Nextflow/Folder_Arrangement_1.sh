#!/bin/bash

# Define the base directory as the current working directory
base_dir="$(pwd)"

# Define the working directory
work_dir="work"

# Function to move directories if they exist
move_directory() {
    local dir_name="$1"
    local source_dir

    # Find the directory
    source_dir=$(find "$work_dir" -type d -name "$dir_name" -print -quit)

    if [ -n "$source_dir" ]; then
        echo "Moving directory: $source_dir to $base_dir"
        mv "$source_dir" "$base_dir/"
    else
        echo "Directory $dir_name not found in $work_dir"
    fi
}

# Move StrelkaSomaticWorkflow directory
move_directory "StrelkaSomaticWorkflow"

# Move StrelkaGermlineWorkflow directory
move_directory "StrelkaGermlineWorkflow"

# Move MantaWorkflow directory
move_directory "MantaWorkflow"

echo "Directories have been moved to the base directory."

