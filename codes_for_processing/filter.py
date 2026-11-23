from pathlib import Path

import pandas as pd


def prompt_path(prompt):
    value = input(f"{prompt}: ").strip()
    if not value:
        raise ValueError(f"No path provided for: {prompt}")
    return Path(value)


try:
    input_path = prompt_path("Enter input TSV path")
    filtered_output_path = prompt_path("Enter filtered TSV output path")
    accession_output_path = prompt_path("Enter assembly accession output txt path")
except ValueError as e:
    print(f"Error: {e}")
    raise SystemExit(1)

if not input_path.exists():
    print(f"Error: input file not found: {input_path}")
    raise SystemExit(1)

filtered_output_path.parent.mkdir(parents=True, exist_ok=True)
accession_output_path.parent.mkdir(parents=True, exist_ok=True)

# 1. Read your NCBI TSV metadata file
df = pd.read_csv(input_path, sep="\t")

# 2. Filter out rows containing ozaenae or rhinoscleromatis
remove_mask = df["Organism Name"].str.contains(
    "ozaenae|rhinoscleromatis", case=False, na=False
)
removed_df = df.loc[remove_mask, ["Organism Name", "Assembly Accession"]]
cleaned_df = df.loc[~remove_mask]

# 3. Save the cleaned results to a new TSV file
cleaned_df.to_csv(filtered_output_path, sep="\t", index=False)

# 4. Save only the assembly accessions to a text file (one per line)
accessions = removed_df["Assembly Accession"].dropna().astype(str)
accessions.to_csv(accession_output_path, index=False, header=False)

print(f"Original entries: {len(df)}")
print(f"Removed entries: {len(removed_df)}")
print("Removed organisms and their Assembly Accession:")
if removed_df.empty:
    print("None")
else:
    print(removed_df.to_string(index=False))
print(f"Cleaned entries kept: {len(cleaned_df)}")
print(f"Assembly accessions saved to {accession_output_path}")