import os
from collections import defaultdict
from Bio import SeqIO

def prompt_text(message):
    while True:
        value = input(message).strip()
        if value:
            return value
        print("A value is required.")


def prompt_float(message):
    while True:
        value = prompt_text(message)
        try:
            return float(value)
        except ValueError:
            print("Please enter a valid number.")


def main():
    blast_tsv = prompt_text("Path to the BLAST TSV file: ")
    reference_fasta = prompt_text("Path to the reference protein FASTA: ")
    query_fasta = prompt_text("Path to the query FASTA: ")
    output_dir = prompt_text("Directory for summary output files: ")
    query_dir = prompt_text("Directory for individual reference FASTA files: ")
    min_pident = prompt_float("Minimum percent identity: ")
    min_qcovs = prompt_float("Minimum query coverage percentage: ")
    max_evalue = prompt_float("Maximum E-value: ")

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(query_dir, exist_ok=True)

    reference_sequences = SeqIO.index(reference_fasta, "fasta")
    print(f"Indexed {len(reference_sequences)} reference sequences.")

    query_hits = defaultdict(list)
    with open(blast_tsv, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            cols = line.split("\t")
            if len(cols) < 13:
                continue

            qseqid = cols[0]
            uniprot_id = qseqid.split("|")[1] if "|" in qseqid else qseqid
            sseqid = cols[1]

            try:
                pident, evalue = float(cols[2]), float(cols[10])
                bitscore, qcovs = float(cols[11]), float(cols[12])
            except ValueError:
                continue

            if pident >= min_pident and qcovs >= min_qcovs and evalue <= max_evalue:
                query_hits[uniprot_id].append((sseqid, pident, qcovs, evalue, bitscore))

    no_reference_hits = []
    mapping = []
    written_ref_files = set()
    selected_ref_ids = set()
    all_queries = [
        record.id.split("|")[1] if "|" in record.id else record.id
        for record in SeqIO.parse(query_fasta, "fasta")
    ]

    for uniprot_id in all_queries:
        hits = query_hits.get(uniprot_id, [])
        if not hits:
            no_reference_hits.append(uniprot_id)
            continue

        hits.sort(key=lambda x: (-x[2], -x[1], -x[4], x[3]))
        best_hit = hits[0]
        reference_id = best_hit[0]

        if reference_id not in reference_sequences:
            no_reference_hits.append(uniprot_id)
            continue

        mapping.append((uniprot_id, reference_id, *best_hit[1:]))
        selected_ref_ids.add(reference_id)

        safe_ref_id = "".join(c if c.isalnum() or c in "._-" else "_" for c in reference_id)
        output_fasta = os.path.join(query_dir, f"ref_{safe_ref_id}.fasta")
        if output_fasta not in written_ref_files:
            SeqIO.write(reference_sequences[reference_id], output_fasta, "fasta")
            written_ref_files.add(output_fasta)

    with open(os.path.join(output_dir, "no_reference_hit.txt"), "w") as f:
        f.write("\n".join(no_reference_hits))

    with open(os.path.join(output_dir, "reference_map.tsv"), "w") as f:
        f.write("query_uniprot\treference_id\tpident\tqcovs\tevalue\tbitscore\n")
        for m in mapping:
            f.write(f"{m[0]}\t{m[1]}\t{m[2]:.2f}\t{m[3]:.2f}\t{m[4]:.2e}\t{m[5]:.1f}\n")

    with open(os.path.join(output_dir, "selected_references.fasta"), "w") as f:
        for rid in sorted(selected_ref_ids):
            SeqIO.write(reference_sequences[rid], f, "fasta")

    print(f"Mapped {len(mapping)} queries. {len(no_reference_hits)} had no hits.")
    print(f"Single-sequence query FASTAs saved in: {query_dir}")


if __name__ == "__main__":
    main()