import math
from pathlib import Path

import pandas as pd

input_dir = Path(input("Enter the folder containing mutation CSV files: ").strip())
summary_csv = Path(input("Enter the mutation summary CSV path: ").strip())
output_dir = Path(input("Enter the output directory path: ").strip())
output_dir.mkdir(parents=True, exist_ok=True)

mutation_files = sorted(input_dir.glob("*_mutations.csv"))
if not mutation_files:
    raise FileNotFoundError(f"No *_mutations.csv files found in: {input_dir}")
if not summary_csv.is_file():
    raise FileNotFoundError(f"Mutation summary CSV not found: {summary_csv}")

summary_df = pd.read_csv(summary_csv)
required_summary_columns = {"MSA File", "Number of Mutations", "Substitutions"}
missing_columns = required_summary_columns - set(summary_df.columns)
if missing_columns:
    raise ValueError(
        "Mutation summary CSV is missing columns: "
        + ", ".join(sorted(missing_columns))
    )

summary_df["summary_key"] = (
    summary_df["MSA File"].astype(str).map(lambda value: Path(value).stem)
)
if summary_df["summary_key"].duplicated().any():
    duplicate_keys = sorted(
        summary_df.loc[
            summary_df["summary_key"].duplicated(keep=False), "summary_key"
        ].unique()
    )
    raise ValueError(
        "Mutation summary CSV contains duplicate MSA File entries: "
        + ", ".join(duplicate_keys)
    )

for input_csv in mutation_files:
    df = pd.read_csv(input_csv)
    required_mutation_columns = {"Mutation Type", "Aligned Position", "Count"}
    missing_columns = required_mutation_columns - set(df.columns)
    if missing_columns:
        raise ValueError(
            f"{input_csv} is missing columns: "
            + ", ".join(sorted(missing_columns))
        )

    mutation_key = input_csv.name.removesuffix("_mutations.csv")
    matching_summary = summary_df[summary_df["summary_key"] == mutation_key]
    if matching_summary.empty:
        raise ValueError(
            f"No matching MSA File found in {summary_csv} for {input_csv.name}"
        )

    summary_row = matching_summary.iloc[0]
    total_mutations = pd.to_numeric(
        summary_row["Number of Mutations"], errors="coerce"
    )
    total_substitutions = pd.to_numeric(
        summary_row["Substitutions"], errors="coerce"
    )
    if (
        pd.isna(total_mutations)
        or pd.isna(total_substitutions)
        or total_mutations <= 0
        or total_substitutions <= 0
        or total_substitutions > total_mutations
    ):
        raise ValueError(
            f"Invalid mutation totals for {input_csv.name}: "
            f"total mutations={summary_row['Number of Mutations']}, "
            f"substitutions={summary_row['Substitutions']}"
        )

    subs_df = df[df["Mutation Type"] == "Substitution"].copy()
    if subs_df.empty:
        raise ValueError(
            f"No substitution rows found in mutation file: {input_csv}"
        )

    subs_df["Count"] = pd.to_numeric(subs_df["Count"], errors="coerce")
    if subs_df["Count"].isna().any() or (subs_df["Count"] < 0).any():
        raise ValueError(f"Invalid substitution counts found in: {input_csv}")

    position_totals = subs_df.groupby("Aligned Position")["Count"].transform("sum")
    subs_df["PSS"] = subs_df["Count"] / position_totals * 100

    global_substitution_frequency = total_substitutions / total_mutations
    subs_df["Log_PSS"] = subs_df["Count"].div(position_totals).map(
        lambda value: math.log(value / global_substitution_frequency)
    )

    output_name = input_csv.name.replace("_mutations.csv", "_PSSM.csv")
    output_csv = output_dir / output_name
    subs_df.to_csv(output_csv, index=False)
    print(f"Substitution data with PSS and Log_PSS saved to: {output_csv}")