import argparse
import csv
from pathlib import Path


def count_hits(csv_folder: Path) -> list[tuple[str, int]]:
    """Return the number of filtered hits in each query CSV file."""
    counts = []

    for csv_path in sorted(csv_folder.glob("*_hit_defs.csv")):
        with csv_path.open(newline="") as csv_file:
            reader = csv.DictReader(csv_file)
            hit_count = sum(1 for row in reader if row.get("Hit_Def"))

        query_name = csv_path.name.removesuffix("_hit_defs.csv")
        counts.append((query_name, hit_count))

    return counts


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Count filtered BLAST hits for each query CSV file."
    )
    parser.add_argument(
        "csv_folder",
        nargs="?",
        type=Path,
        default=Path("results/hit_filtered_results/filtered_hits_csv"),
        help="Folder containing *_hit_defs.csv files.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Optional path for a CSV summary file.",
    )
    args = parser.parse_args()

    if not args.csv_folder.is_dir():
        parser.error(f"CSV folder not found: {args.csv_folder}")

    counts = count_hits(args.csv_folder)
    if not counts:
        print(f"No *_hit_defs.csv files found in: {args.csv_folder}")
        return

    print("Query\tFiltered_Hits")
    for query_name, hit_count in counts:
        print(f"{query_name}\t{hit_count}")

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", newline="") as output_file:
            writer = csv.writer(output_file)
            writer.writerow(["Query", "Filtered_Hits"])
            writer.writerows(counts)
        print(f"\nSummary saved to: {args.output}")


if __name__ == "__main__":
    main()