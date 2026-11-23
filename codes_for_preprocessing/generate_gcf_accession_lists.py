#!/usr/bin/env python3
"""Generate cleaned and skipped GCF accession lists from NCBI dataset TSV files."""

import argparse
import csv
import sys
from pathlib import Path


ACCESSION_COLUMN = "Assembly Accession"


def read_accessions(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8", newline="") as input_file:
        reader = csv.DictReader(input_file, delimiter="\t")
        if not reader.fieldnames or ACCESSION_COLUMN not in reader.fieldnames:
            raise ValueError(f"{path} must contain a '{ACCESSION_COLUMN}' column")

        accessions = []
        seen = set()
        for row in reader:
            accession = (row.get(ACCESSION_COLUMN) or "").strip()
            if accession and accession not in seen:
                seen.add(accession)
                accessions.append(accession)
        return accessions


def write_accessions(path: Path, accessions: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as output_file:
        output_file.write("\n".join(accessions))
        if accessions:
            output_file.write("\n")


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Create a list of cleaned GCF accessions and a list of accessions "
            "to skip because they are absent from the cleaned TSV."
        )
    )
    parser.add_argument("ncbi_dataset", type=Path, help="Original ncbi_dataset.tsv")
    parser.add_argument("cleaned_dataset", type=Path, help="Filtered cleaned_ncbi_dataset.tsv")
    parser.add_argument(
        "--cleaned-output",
        type=Path,
        default=Path("cleaned_gcf_accessions.txt"),
        help="Output file for accessions present in the cleaned TSV",
    )
    parser.add_argument(
        "--skip-output",
        type=Path,
        default=Path("skipped_gcf_accessions.txt"),
        help="Output file for accessions absent from the cleaned TSV",
    )
    args = parser.parse_args()

    for path in (args.ncbi_dataset, args.cleaned_dataset):
        if not path.is_file():
            print(f"Input file not found: {path}", file=sys.stderr)
            return 1

    try:
        all_accessions = read_accessions(args.ncbi_dataset)
        cleaned_accessions = read_accessions(args.cleaned_dataset)
    except (OSError, UnicodeError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1

    cleaned_set = set(cleaned_accessions)
    skipped_accessions = [
        accession for accession in all_accessions if accession not in cleaned_set
    ]

    try:
        write_accessions(args.cleaned_output, cleaned_accessions)
        write_accessions(args.skip_output, skipped_accessions)
    except OSError as error:
        print(f"Error writing output: {error}", file=sys.stderr)
        return 1

    print(f"Wrote {args.cleaned_output} ({len(cleaned_accessions)} accessions)")
    print(f"Wrote {args.skip_output} ({len(skipped_accessions)} accessions to skip)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
