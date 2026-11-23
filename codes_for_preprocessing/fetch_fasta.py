import requests
from pathlib import Path
import time
import sys


def prompt_text(prompt: str) -> str:
    return input(prompt + " ").strip()


def read_ids(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8") as fh:
        seen = set()
        ids = []
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            if s not in seen:
                seen.add(s)
                ids.append(s)
    return ids


def extract_first_fasta(fasta_text: str) -> str | None:
    """Extracts only the first complete FASTA record if multiple entries are returned."""
    stripped = fasta_text.strip()
    if not stripped.startswith(">"):
        return None

    # Split records by the start of subsequent entries
    records = stripped.split("\n>")
    first_record = records[0]
    if not first_record.startswith(">"):
        first_record = ">" + first_record

    return first_record + "\n"


def fetch_with_retries(
    session: requests.Session, uid: str, tax_id: str = "", attempts: int = 3
) -> str | None:
    # Build Tier 2 URL dynamically based on whether taxonomy filtering is provided
    if tax_id:
        tier2_url = f"https://rest.uniprot.org/uniprotkb/stream?query=sec_acc:{uid}+AND+taxonomy_id:{tax_id}&format=fasta"
    else:
        tier2_url = f"https://rest.uniprot.org/uniprotkb/stream?query=sec_acc:{uid}&format=fasta"

    # 3-tier fallback sequence
    fallback_urls = [
        # Tier 1: Direct active UniProtKB entry
        f"https://rest.uniprot.org/uniprotkb/{uid}.fasta",
        # Tier 2: Secondary / merged UniProtKB entry
        tier2_url,
        # Tier 3: UniParc sequence archive
        f"https://rest.uniprot.org/uniparc/stream?query={uid}&format=fasta",
    ]

    for url in fallback_urls:
        backoff = 1.0
        for attempt in range(1, attempts + 1):
            try:
                resp = session.get(url, timeout=30)
            except Exception as exc:
                print(f"[{attempt}/{attempts}] Connection error on {uid}: {exc}")
                if attempt < attempts:
                    time.sleep(backoff)
                    backoff *= 2
                    continue
                break

            # If an entry is found and contains text
            if resp.status_code == 200 and resp.text.strip():
                single_record = extract_first_fasta(resp.text)
                if single_record:
                    return single_record

            # Rate limits (429) or server errors (5xx) -> retry current endpoint
            if resp.status_code == 429 or 500 <= resp.status_code < 600:
                retry_after = resp.headers.get("Retry-After")
                wait = int(retry_after) if retry_after and retry_after.isdigit() else backoff
                print(f"[{attempt}/{attempts}] Server busy ({resp.status_code}) for {uid}. Retrying...")
                time.sleep(wait)
                backoff *= 2
                continue

            # 404 Not Found or empty search result -> proceed to the next fallback tier
            break

    return None


def main() -> None:
    inp = prompt_text("Path to UniProt IDs txt file:")
    outp = prompt_text("Output FASTA path to write combined sequences:")
    tax_filter = prompt_text(
        "Optional NCBI Taxonomy ID to filter (e.g., 573, or press Enter to skip):"
    )

    inpath = Path(inp)
    outpath = Path(outp)
    if not inpath.exists():
        print(f"Input file does not exist: {inpath}", file=sys.stderr)
        sys.exit(1)

    outpath.parent.mkdir(parents=True, exist_ok=True)

    ids = read_ids(inpath)
    if not ids:
        print("No UniProt IDs found in input.")
        return

    session = requests.Session()
    failed_path = outpath.parent / "failed_uniprot_ids.txt"

    with outpath.open("w", encoding="utf-8") as outfh, failed_path.open("w", encoding="utf-8") as failfh:
        for i, uid in enumerate(ids, start=1):
            fasta_entry = fetch_with_retries(session, uid, tax_id=tax_filter, attempts=3)
            if fasta_entry:
                # Prepend the queried accession to the FASTA header
                lines = fasta_entry.splitlines()
                original_header = lines[0][1:]  # strip leading '>'
                lines[0] = f">{uid} | {original_header}"
                standardized_entry = "\n".join(lines) + "\n"

                outfh.write(standardized_entry)
                print(f"[{i}/{len(ids)}] Downloaded: {uid}")
            else:
                failfh.write(uid + "\n")
                print(f"[{i}/{len(ids)}] Failed across all tiers: {uid}")

            # Polite delay to prevent rate limiting
            time.sleep(0.2)


if __name__ == "__main__":
    main()