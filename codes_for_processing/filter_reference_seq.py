import csv
import json

MIN_IDENTITY = 0.90  # 90%
MIN_COVERAGE = 0.70  # 70%
MAX_EVALUE = 1e-5    # 1e-5

raw_blast_fmt6 = "results/pdb_to_ref_raw.fmt6"
output_json = "results/ref_map.json"
output_csv = "results/ref_map_details.csv"

best_hits = {}

with open(raw_blast_fmt6, "r") as f:
    reader = csv.reader(f, delimiter="\t")
    for row in reader:
        if not row:
            continue
            
        qseqid = row[0]
        sseqid = row[1]
        pident = float(row[2]) / 100.0  # Convert to fraction
        align_len = float(row[3])
        qlen = float(row[4])
        nident = int(row[5])
        evalue = float(row[6])

        # Calculate coverage relative to query length
        coverage = align_len / qlen

        # Apply filtering threshold
        if pident >= MIN_IDENTITY and coverage >= MIN_COVERAGE and evalue < MAX_EVALUE:
            if qseqid not in best_hits or evalue < best_hits[qseqid]["evalue"]:
                best_hits[qseqid] = {
                    "ref_id": sseqid,
                    "pident": pident,
                    "coverage": coverage,
                    "evalue": evalue,
                    "align_len": align_len,
                    "qlen": qlen,
                    "nident": nident
                }

# Create 1-to-1 query -> reference mapping dictionary
ref_mapping = {qseqid: hit_data["ref_id"] for qseqid, hit_data in best_hits.items()}

with open(output_json, "w") as jout:
    json.dump(ref_mapping, jout, indent=2)
with open(output_csv, "w", newline='') as cout:
    writer = csv.writer(cout)
    writer.writerow(["qseqid", "ref_id", "pident_percent", "coverage", "qlen", "align_len", "nident", "evalue"])
    for qseqid, hit in best_hits.items():
        writer.writerow([
            qseqid,
            hit["ref_id"],
            round(hit["pident"] * 100, 6),
            hit["coverage"],
            int(hit["qlen"]),
            int(hit["align_len"]),
            int(hit["nident"]),
            hit["evalue"]
        ])

print(f"✅ Created {output_json} and {output_csv} with {len(ref_mapping)} mapped reference sequences!")