#!/usr/bin/env python3
"""Keep Klebsiella pneumoniae records except ozaenae and rhinoscleromatis subspecies."""

import argparse
import csv
import sys
from pathlib import Path


ORGANISM_COLUMN = "Organism Name"
EXCLUDED_SUBSPECIES = (
    "subsp. ozaenae",
    "subsp. rhinoscleromatis",
)


def is_included_organism(organism: str) -> bool:
    """Return whether an organism should remain in the cleaned dataset."""
    normalized = organism.strip().casefold()
    return (
        normalized.startswith("klebsiella pneumoniae")
        and not any(subspecies in normalized for subspecies in EXCLUDED_SUBSPECIES)
    )


def filter_dataset(input_path: Path, output_path: Path) -> tuple[int, int]:
    with input_path.open("r", encoding="utf-8", newline="") as input_file:
        reader = csv.DictReader(input_file, delimiter="\t")
        if not reader.fieldnames or ORGANISM_COLUMN not in reader.fieldnames:
            raise ValueError(f"{input_path} must contain an '{ORGANISM_COLUMN}' column")

        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("w", encoding="utf-8", newline="") as output_file:
            writer = csv.DictWriter(
                output_file,
                fieldnames=reader.fieldnames,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()

            total = included = 0
            for row in reader:
                total += 1
                if is_included_organism(row.get(ORGANISM_COLUMN, "")):
                    writer.writerow(row)
                    included += 1

    return total, included


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create cleaned_ncbi_dataset.tsv by filtering organism names."
    )
    parser.add_argument("input_tsv", type=Path, help="Original ncbi_dataset.tsv")
    parser.add_argument(
        "output_tsv",
        type=Path,
        help="Output cleaned TSV, for example cleaned_ncbi_dataset.tsv",
    )
    args = parser.parse_args()

    if not args.input_tsv.is_file():
        print(f"Input file not found: {args.input_tsv}", file=sys.stderr)
        return 1

    try:
        total, included = filter_dataset(args.input_tsv, args.output_tsv)
    except (OSError, UnicodeError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1

    print(
        f"Wrote {args.output_tsv}: kept {included} of {total} records "
        f"({total - included} excluded)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
