from Bio import SeqIO
import pandas as pd
import os


# Function to determine conservative/non-conservative substitution
def get_substitution_type(original, mutated):
    conservative_groups = [
        {'A', 'V', 'L', 'I', 'M'},      # Aliphatic
        {'F', 'Y', 'W'},                # Aromatic
        {'D', 'E'},                     # Acidic
        {'K', 'R', 'H'},                # Basic
        {'S', 'T', 'N', 'Q'},            # Polar uncharged
        {'G', 'P', 'C'}                 # Special cases
    ]

    for group in conservative_groups:
        if original in group and mutated in group:
            return "Conservative"

    return "Non-Conservative"


# Map alignment positions to positions in the original ungapped
# reference sequence
def map_alignment_to_reference_positions(reference_sequence):
    position_map = []
    reference_pos = 0

    for aa in reference_sequence:
        if aa == "-":
            position_map.append(None)
        else:
            reference_pos += 1
            position_map.append(reference_pos)

    return position_map


# Main processing function
def extract_mutations(msa_folder, output_root_folder):

    mutation_folder = os.path.join(output_root_folder, "mutations")
    os.makedirs(mutation_folder, exist_ok=True)

    summary_data = []

    # Process every MSA file
    for msa_filename in os.listdir(msa_folder):

        if not msa_filename.endswith(".faa"):
            continue

        msa_path = os.path.join(msa_folder, msa_filename)

        # ---------------------------------------------------------
        # Read MSA
        # ---------------------------------------------------------
        with open(msa_path, "r") as handle:
            records = list(SeqIO.parse(handle, "fasta"))

        if not records:
            print(f"Skipped {msa_filename}: empty file")
            continue

        # ---------------------------------------------------------
        # FIRST SEQUENCE = REFERENCE SEQUENCE
        # ---------------------------------------------------------
        reference_sequence = str(records[0].seq)

        aligned_sequences = [
            str(record.seq)
            for record in records
        ]

        # ---------------------------------------------------------
        # Check that all sequences have the same alignment length
        # ---------------------------------------------------------
        alignment_lengths = {
            len(seq) for seq in aligned_sequences
        }

        if len(alignment_lengths) != 1:
            print(
                f"Skipped {msa_filename}: "
                f"sequences have different alignment lengths"
            )
            continue

        # ---------------------------------------------------------
        # Total homologs
        # ---------------------------------------------------------
        total_sequences = len(aligned_sequences)
        total_homologs = total_sequences - 1

        if total_homologs <= 0:
            print(
                f"Skipped {msa_filename}: "
                f"no homolog sequences found"
            )
            continue

        # ---------------------------------------------------------
        # Map alignment columns to reference positions
        # ---------------------------------------------------------
        reference_positions = map_alignment_to_reference_positions(
            reference_sequence
        )

        # ---------------------------------------------------------
        # Count amino acids ONLY in homolog sequences
        # ---------------------------------------------------------
        aa_matrix = [
            {}
            for _ in range(len(reference_sequence))
        ]

        for sequence in aligned_sequences[1:]:

            for index, aa in enumerate(sequence):

                aa_matrix[index][aa] = (
                    aa_matrix[index].get(aa, 0) + 1
                )

        # ---------------------------------------------------------
        # Store mutation information here
        # ---------------------------------------------------------
        mutations = []

        # ---------------------------------------------------------
        # Examine every alignment position
        # ---------------------------------------------------------
        for index, count_dict in enumerate(aa_matrix):

            reference_residue = reference_sequence[index]

            aligned_position = index + 1
            reference_position = reference_positions[index]

            # =====================================================
            # CASE 1: REFERENCE HAS AN AMINO ACID
            # =====================================================
            if reference_residue != "-":

                # -------------------------------------------------
                # DELETION
                # Reference = amino acid
                # Homolog   = "-"
                # -------------------------------------------------
                if "-" in count_dict:

                    deletion_count = count_dict["-"]

                    percentage = round(
                        (deletion_count / total_homologs) * 100,
                        2
                    )

                    mutations.append({
                        "Aligned Position": aligned_position,
                        "Reference Position": reference_position,
                        "Reference Residue": reference_residue,
                        "Mutation Type": "Deletion",
                        "Mutated Residue": "",
                        "Count": deletion_count,
                        "Total Homologs": total_homologs,
                        "Percentage": percentage,
                        "Substitution Type": ""
                    })

                # -------------------------------------------------
                # SUBSTITUTION
                # Reference = amino acid
                # Homolog   = different amino acid
                # -------------------------------------------------
                for mutated_residue, count in count_dict.items():

                    # Ignore gaps and same residue
                    if mutated_residue == "-":
                        continue

                    if mutated_residue == reference_residue:
                        continue

                    percentage = round(
                        (count / total_homologs) * 100,
                        2
                    )

                    substitution_type = get_substitution_type(
                        reference_residue,
                        mutated_residue
                    )

                    mutations.append({
                        "Aligned Position": aligned_position,
                        "Reference Position": reference_position,
                        "Reference Residue": reference_residue,
                        "Mutation Type": "Substitution",
                        "Mutated Residue": mutated_residue,
                        "Count": count,
                        "Total Homologs": total_homologs,
                        "Percentage": percentage,
                        "Substitution Type": substitution_type
                    })

            # =====================================================
            # CASE 2: REFERENCE HAS A GAP
            # =====================================================
            else:

                # -------------------------------------------------
                # INSERTION
                # Reference = "-"
                # Homolog   = amino acid
                # -------------------------------------------------
                for inserted_residue, count in count_dict.items():

                    if inserted_residue == "-":
                        continue

                    percentage = round(
                        (count / total_homologs) * 100,
                        2
                    )

                    mutations.append({
                        "Aligned Position": aligned_position,
                        "Reference Position": "N/A",
                        "Reference Residue": "-",
                        "Mutation Type": "Insertion",
                        "Mutated Residue": inserted_residue,
                        "Count": count,
                        "Total Homologs": total_homologs,
                        "Percentage": percentage,
                        "Substitution Type": ""
                    })

        # ---------------------------------------------------------
        # Convert mutation list to Pandas DataFrame
        # ---------------------------------------------------------
        mutation_df = pd.DataFrame(
            mutations,
            columns=[
                "Aligned Position",
                "Reference Position",
                "Reference Residue",
                "Mutation Type",
                "Mutated Residue",
                "Count",
                "Total Homologs",
                "Percentage",
                "Substitution Type"
            ]
        )

        # ---------------------------------------------------------
        # Save mutation DataFrame
        # ---------------------------------------------------------
        if not mutation_df.empty:

            mutation_counts = mutation_df["Mutation Type"].value_counts()
            deletion_count = int(mutation_counts.get("Deletion", 0))
            insertion_count = int(mutation_counts.get("Insertion", 0))
            substitution_count = int(mutation_counts.get("Substitution", 0))

            output_filename = msa_filename.replace(
                ".faa",
                "_mutations.csv"
            )

            output_path = os.path.join(
                mutation_folder,
                output_filename
            )

            mutation_df.to_csv(
                output_path,
                index=False
            )

            summary_data.append({
                "MSA File": msa_filename,
                "Status": "Mutations found",
                "Number of Mutations": len(mutation_df),
                "Deletions": deletion_count,
                "Insertions": insertion_count,
                "Substitutions": substitution_count,
                "Total Homologs": total_homologs
            })

            print(
                f"Processed: {msa_filename} → "
                f"{output_filename}"
            )

        else:

            summary_data.append({
                "MSA File": msa_filename,
                "Status": "No mutations found",
                "Number of Mutations": 0,
                "Deletions": 0,
                "Insertions": 0,
                "Substitutions": 0,
                "Total Homologs": total_homologs
            })

            print(
                f"Processed: {msa_filename} → "
                f"no mutations"
            )

    # -------------------------------------------------------------
    # Create summary DataFrame
    # -------------------------------------------------------------
    summary_df = pd.DataFrame(
        summary_data,
        columns=[
            "MSA File",
            "Status",
            "Number of Mutations",
            "Deletions",
            "Insertions",
            "Substitutions",
            "Total Homologs"
        ]
    )

    # -------------------------------------------------------------
    # Save summary
    # -------------------------------------------------------------
    summary_path = os.path.join(
        output_root_folder,
        "mutation_summary.csv"
    )

    summary_df.to_csv(
        summary_path,
        index=False
    )

    print(
        f"\nMutation summary written to: {summary_path}"
    )


# ================================================================
# Run script
# ================================================================
if __name__ == "__main__":

    msa_input_folder = input(
        "Enter path to your folder containing MSA (.faa) files: "
    ).strip()

    output_folder = input(
        "Enter path for output folder to store mutation results: "
    ).strip()

    extract_mutations(
        msa_input_folder,
        output_folder
    )