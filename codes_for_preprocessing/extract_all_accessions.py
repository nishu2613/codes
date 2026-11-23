#!/usr/bin/env python3
"""Extract all unique UniProt accessions from a CIF accession CSV."""
from pathlib import Path
import csv
import sys


def main() -> None:
    input_text = input("Input CIF accession CSV: ").strip()
    output_text = input("Output TXT (one UniProt accession per line): ").strip()

    if not input_text or not output_text:
        print("Both input and output paths are required.", file=sys.stderr)
        sys.exit(1)

    input_path = Path(input_text)
    output_path = Path(output_text)

    if not input_path.is_file():
        print(f"Input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    accessions = []
    seen = set()

    with input_path.open("r", newline="", encoding="utf-8") as input_file:
        reader = csv.DictReader(input_file)
        if not reader.fieldnames or "accessions" not in reader.fieldnames:
            print("Input CSV must contain an 'accessions' column.", file=sys.stderr)
            sys.exit(1)

        for row in reader:
            value = (row.get("accessions") or "").strip()
            if not value:
                continue

            for accession in value.split(";"):
                accession = accession.strip()
                if accession and accession not in {"NIL", ".", "?"} and accession not in seen:
                    seen.add(accession)
                    accessions.append(accession)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as output_file:
        for accession in accessions:
            output_file.write(f"{accession}\n")

    print(f"Wrote {output_path} ({len(accessions)} unique accessions)")


if __name__ == "__main__":
    main()
