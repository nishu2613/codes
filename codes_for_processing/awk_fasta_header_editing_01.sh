#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <path_to_folder_containing_GCF_directories> <path_to_skip_list_txt>"
    echo "Example: $0 /home/neelamnishanth/Kp_data_processing/Genomes_data/genomes/ncbi_dataset/data /home/neelamnishanth/skip_gcf.txt"
    exit 1
fi

input_path="$1"
skip_list_file="$2"
output_file="all_genomes.faa"

if [ ! -d "$input_path" ]; then
    echo "Error: directory not found: $input_path"
    exit 1
fi

if [ ! -f "$skip_list_file" ]; then
    echo "Error: skip list file not found: $skip_list_file"
    exit 1
fi

if ! compgen -G "$input_path"/GCF* >/dev/null; then
    echo "Error: no GCF* folders found inside: $input_path"
    exit 1
fi

# Build skip map from the user-provided txt file
declare -A skip_map=()
while IFS= read -r line || [ -n "$line" ]; do
    clean_line="${line%$'\r'}"
    clean_line="${clean_line//[$' \t\r\n']/}"
    if [ -n "$clean_line" ]; then
        skip_map["$clean_line"]=1
    fi
done < "$skip_list_file"

# Build final list of directories to process
all_dirs=("$input_path"/GCF*)
filtered_dirs=()
for dir in "${all_dirs[@]}"; do
    [ -d "$dir" ] || continue

    gcf_id=$(basename "$dir")
    if [ -n "${skip_map[$gcf_id]:-}" ]; then
        echo "Skipping GCF folder: $gcf_id"
        continue
    fi

    filtered_dirs+=("$dir")
done

if [ "${#filtered_dirs[@]}" -eq 0 ]; then
    echo "Error: no eligible GCF folders remain after applying the skip list."
    exit 1
fi

> "$output_file"

for dir in "${filtered_dirs[@]}"; do
    gcf_id=$(basename "$dir")

    for file in "$dir"/*protein.faa; do
        [ -f "$file" ] || continue

        awk -v gcf="$gcf_id" '
        /^>/ {
            line = substr($0, 2)
            split(line, arr, " ")
            protein_acc = arr[1]
            rest = substr(line, length(protein_acc) + 2)
            print ">" protein_acc "_" gcf " " gcf " " protein_acc " " rest
            next
        }
        { print }
        ' "$file"
    done
done >> "$output_file"

echo "Merged FASTA saved as: $output_file"
