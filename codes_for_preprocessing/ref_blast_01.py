import os
import subprocess


def prompt_required(message: str) -> str:
    while True:
        value = input(message).strip()
        if value:
            return value
        print("A value is required.")


def main() -> None:
    input_fasta = prompt_required("Path to the input FASTA file: ")
    blast_db = prompt_required("Path to the BLAST reference database: ")
    output_dir = prompt_required("Directory for the BLAST output: ")
    output_filename = prompt_required("Output TSV filename: ")
    evalue = prompt_required("E-value threshold: ")
    max_target_seqs = prompt_required("Maximum target sequences per query: ")
    num_threads = prompt_required("Number of BLAST threads: ")

    os.makedirs(output_dir, exist_ok=True)
    output_tsv = os.path.join(output_dir, output_filename)

    command = [
        "blastp",
        "-query", input_fasta,
        "-db", blast_db,
        "-outfmt", "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qcovs",
        "-evalue", evalue,
        "-max_target_seqs", max_target_seqs,
        "-num_threads", num_threads,
        "-out", output_tsv,
    ]

    print("Starting single-batch BLAST search... This may take a moment.")
    result = subprocess.run(command, capture_output=True, text=True)

    if result.returncode != 0:
        print("ERROR: BLAST failed")
        print(result.stderr)
    else:
        print(f"BLAST complete! All results saved to: {output_tsv}")


if __name__ == "__main__":
    main()