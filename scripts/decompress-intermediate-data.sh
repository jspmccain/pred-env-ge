#!/bin/bash

# Define the folder containing the compressed files
input_folder="../data/tester_int/" # Adjust this path if needed

# Loop through all the .7z.001 files in the folder
for header_file in "$input_folder"/*.7z.001; do
    # Skip if no files match the pattern
    if [[ ! -f "$header_file" ]]; then
        echo "No .7z.001 files found in the folder."
        break
    fi

    # Extract the base name without the extension
    base_name=$(basename -- "$header_file" .001)
    
    # Decompress and merge the parts
    7z x "$header_file"
    
    # Notify the user
    echo "Decompressed and merged parts for $base_name"
done

# Move the rds files out of /scripts and into intermediate_data
mv *rds $input_folder

