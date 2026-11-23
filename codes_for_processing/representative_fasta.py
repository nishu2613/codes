import os

fasta_folder = "structure_data/kp_fasta_files"  # Folder containing individual .fasta files
rep_list_file = "results/representatives.txt"
output_combined_fasta = "results/representatives.fasta"

with open(rep_list_file, "r") as f:
    rep_files = [line.strip() for line in f if line.strip()]

with open(output_combined_fasta, "w") as outfile:
    for file_name in rep_files:
        file_path = os.path.join(fasta_folder, file_name)
        if os.path.exists(file_path):
            with open(file_path, "r") as infile:
                outfile.write(infile.read().strip() + "\n")

print(f"✅ Created {output_combined_fasta} with {len(rep_files)} representative sequences!")