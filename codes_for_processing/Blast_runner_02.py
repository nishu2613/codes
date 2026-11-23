import os
import subprocess

# -------- USER INPUTS --------
fasta_folder = input("Enter path to the FASTA query folder: ").strip()
blast_db = input("Enter path to BLAST database: ").strip()
output_folder = input("Enter path to output folder: ").strip()
blast_exec = input("Enter BLAST program (blastp/blastn): ").strip()
threads = input("Enter number of threads to use: ").strip()
max_seqs = input("Enter max target sequences (e.g., 20000 to override the 500-hit limit): ").strip()

# Create output folder if needed
os.makedirs(output_folder, exist_ok=True)
fasta_files = sorted(
    f for f in os.listdir(fasta_folder)
    if f.endswith(".fasta") and os.path.isfile(os.path.join(fasta_folder, f))
)

print(f"\nTotal FASTA entries to process: {len(fasta_files)}\n")

for fasta_name in fasta_files:
    input_path = os.path.join(fasta_folder, fasta_name)
    base_name = os.path.splitext(fasta_name)[0]
    output_path = os.path.join(output_folder, f"{base_name}.xml")
            
    # Skip already processed files (prevents starting over from scratch)
    if os.path.exists(output_path):
        print(f"Skipping {fasta_name} - output already exists.")
        continue

    # Added max_target_seqs and an explicit E-value filter
    cmd = [
        blast_exec,
        "-query", input_path,
        "-db", blast_db,
        "-out", output_path,
        "-outfmt", "5",
        "-evalue", "1e-5",
        "-max_target_seqs", max_seqs,
        "-num_threads", threads
    ]

    print(f"Running BLAST for: {fasta_name}")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"BLAST failed for {fasta_name}: {e}")

print("\nBLAST completed.")