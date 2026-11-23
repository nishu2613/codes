import csv
import json
from collections import defaultdict
from pathlib import Path

import gemmi


def extract_accessions(cif_path: Path) -> list[dict]:
    doc = gemmi.cif.read_file(str(cif_path))
    block = doc.sole_block()

    results = []
    seen = set()

    try:
        table = block.find("_struct_ref.", 
                           ["db_name", "pdbx_db_accession", "entity_id"])
    except RuntimeError:
        return [{"accession": "NIL", "entity_id": "NIL"}]

    for row in table:
        db_name   = row[0].strip()
        accession = row[1].strip()
        entity_id = row[2].strip()

        if db_name != "UNP":
            continue
        if not accession or accession in {".", "?"}:
            continue
            
        # FIX: Track the unique combination of accession AND entity_id
        pair = (accession, entity_id)
        if pair not in seen:
            seen.add(pair)
            results.append({"accession": accession, "entity_id": entity_id})

    if not results:
        return [{"accession": "NIL", "entity_id": "NIL"}]

    return results


def build_clusters(file_to_accessions: dict[str, list[dict]]) -> list[dict]:
	exact_groups = defaultdict(list)
	for filename, accessions_dicts in file_to_accessions.items():
		# Extract just the accession strings for clustering
		accession_strings = [item["accession"] for item in accessions_dicts]
		accession_key = tuple(sorted(set(accession_strings)))
		exact_groups[accession_key].append((filename, accessions_dicts))

	clusters = []
	for accession_key, files_with_dicts in exact_groups.items():
		clusters.append(
			{
				"cluster_id": len(clusters) + 1,
				"file_count": len(files_with_dicts),
				"accession_count": len(accession_key),
				"files": sorted([f[0] for f in files_with_dicts]),
				"accessions": list(accession_key),
				"file_details": [
					{"file": f[0], "accession_data": f[1]}
					for f in files_with_dicts
				],
			}
		)

	clusters.sort(key=lambda c: (-c["file_count"], c["cluster_id"]))
	for index, cluster in enumerate(clusters, start=1):
		cluster["cluster_id"] = index

	return clusters


def write_outputs(
	output_dir: Path,
	file_to_accessions: dict[str, list[dict]],
	clusters: list[dict],
	file_entity_info: dict[str, dict[str, list[str]]],
	classification_file_to_accessions: dict[str, list[dict]],
) -> tuple[Path, Path, Path, Path, Path, Path]:
	output_dir.mkdir(parents=True, exist_ok=True)

	json_path = output_dir / "cif_accession_clusters.json"
	csv_path = output_dir / "cif_file_accessions.csv"
	cluster_csv_path = output_dir / "cif_clusters.csv"
	homopolymer_path = output_dir / "homopolymer_pdb_ids.csv"
	heteropolymer_path = output_dir / "heteropolymer_pdb_ids.csv"
	unmapped_path = output_dir / "unmapped_pdb_ids.csv"

	summary = {
		"total_files": len(file_to_accessions),
		"files_with_accessions": sum(1 for vals in file_to_accessions.values() if vals and vals[0].get("accession") != "NIL"),
		"total_clusters": len(clusters),
		"clusters": clusters,
		"file_accessions": [
			{
				"file": filename,
				"accession_count": len(accessions_dicts),
				"accessions_data": accessions_dicts,
				"accessions": [a.get("accession") for a in accessions_dicts],
				**file_entity_info.get(filename, {}),
			}
			for filename, accessions_dicts in sorted(file_to_accessions.items())
		],
	}

	with json_path.open("w", encoding="utf-8") as handle:
		json.dump(summary, handle, indent=2)

	with csv_path.open("w", newline="", encoding="utf-8") as handle:
		writer = csv.writer(handle)
		writer.writerow([
			"file",
			"total_entity_count",
			"total_entity_ids",
			"accession_entity_count",
			"accession_entity_ids",
			"accession_count",
			"accessions",
			"entity_ids",
		])
		for filename, accessions_dicts in sorted(file_to_accessions.items()):
			accession_strings = [a.get("accession") for a in accessions_dicts]
			entity_ids = [a.get("entity_id") for a in accessions_dicts]
			entity_info = file_entity_info.get(filename, {})
			writer.writerow([
				filename,
				entity_info.get("total_entity_count", ""),
				";".join(entity_info.get("total_entity_ids", [])),
				entity_info.get("accession_entity_count", ""),
				";".join(entity_info.get("accession_entity_ids", [])),
				len(accessions_dicts),
				";".join(accession_strings),
				";".join(entity_ids),
			])

	with cluster_csv_path.open("w", newline="", encoding="utf-8") as handle:
		writer = csv.writer(handle)
		writer.writerow(
			["cluster_id", "file_count", "accession_count", "files", "accessions"]
		)
		for cluster in clusters:
			writer.writerow(
				[
					cluster["cluster_id"],
					cluster["file_count"],
					cluster["accession_count"],
					";".join(cluster["files"]),
					";".join(cluster["accessions"]),
				]
			)

	classified_pdb_ids = {
		"homopolymer": set(),
		"heteropolymer": set(),
		"unmapped": set(),
	}
	classified_accessions = {
		"homopolymer": {},
		"heteropolymer": {},
		"unmapped": {},
	}
	for filename, accessions_dicts in classification_file_to_accessions.items():
		unique_accessions = {
			item.get("accession")
			for item in accessions_dicts
			if item.get("accession") not in {None, "", "NIL", ".", "?"}
		}
		pdb_id = Path(filename).stem.upper()
		if not unique_accessions:
			classified_pdb_ids["unmapped"].add(pdb_id)
			classified_accessions["unmapped"][pdb_id] = set()
		elif len(unique_accessions) == 1:
			classified_pdb_ids["homopolymer"].add(pdb_id)
			classified_accessions["homopolymer"][pdb_id] = unique_accessions
		else:
			classified_pdb_ids["heteropolymer"].add(pdb_id)
			classified_accessions["heteropolymer"][pdb_id] = unique_accessions

	for path, category in (
		(homopolymer_path, "homopolymer"),
		(heteropolymer_path, "heteropolymer"),
		(unmapped_path, "unmapped"),
	):
		with path.open("w", newline="", encoding="utf-8") as handle:
			writer = csv.writer(handle)
			writer.writerow(["pdb_id", "accessions"])
			for pdb_id in sorted(classified_pdb_ids[category]):
				accessions = sorted(classified_accessions[category][pdb_id])
				writer.writerow([pdb_id, ";".join(accessions)])

	return (
		json_path,
		csv_path,
		cluster_csv_path,
		homopolymer_path,
		heteropolymer_path,
		unmapped_path,
	)


def load_allowed_pdb_ids(pdb_id_csv: Path | None) -> set[str]:
	if pdb_id_csv is None:
		return set()

	with pdb_id_csv.open("r", newline="", encoding="utf-8") as handle:
		dreader = csv.DictReader(handle)
		if not dreader.fieldnames or "pdb_id" not in dreader.fieldnames:
			raise ValueError(f"PDB ID CSV must contain a 'pdb_id' column: {pdb_id_csv}")
		return {str(row["pdb_id"]).strip() for row in dreader if str(row.get("pdb_id", "")).strip()}


def load_entity_summary(entity_summary_csv: Path | None) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
	"""Return all entity IDs and biological entity IDs grouped by PDB ID.

	If `entity_summary_csv` is None or the file is missing, returns an empty dict.
	Expects a CSV with columns: `pdb_id`, `entity_id`, `is_biological`.
	"""
	if entity_summary_csv is None:
		return {}, {}

	if not entity_summary_csv.exists():
		return {}

	def _is_true(val: object) -> bool:
		if val is None:
			return False
		s = str(val).strip().lower()
		return s in {"true", "1", "yes", "y", "t"}

	all_entities: dict[str, set[str]] = {}
	biological_entities: dict[str, set[str]] = {}
	with entity_summary_csv.open("r", newline="", encoding="utf-8") as handle:
		dreader = csv.DictReader(handle)
		if not dreader.fieldnames or not {"pdb_id", "entity_id", "is_biological"}.issubset(set(dreader.fieldnames)):
			raise ValueError(f"Entity summary CSV must contain pdb_id, entity_id, and is_biological columns: {entity_summary_csv}")
		for row in dreader:
			pdb_id = str(row.get("pdb_id", "")).strip()
			ent_id = str(row.get("entity_id", "")).strip()
			if not pdb_id or not ent_id:
				continue
			all_entities.setdefault(pdb_id, set()).add(ent_id)
			if _is_true(row.get("is_biological")):
				biological_entities.setdefault(pdb_id, set()).add(ent_id)
	return all_entities, biological_entities


def run(cif_dir: Path, output_dir: Path, pdb_id_csv: Path | None = None, entity_summary_csv: Path | None = None) -> None:
	all_cif_files = sorted(cif_dir.glob("*.cif"))
	if not all_cif_files:
		raise FileNotFoundError(f"No CIF files found in: {cif_dir}")

	allowed_pdb_ids = load_allowed_pdb_ids(pdb_id_csv)
	if allowed_pdb_ids:
		cif_files = [
			cif_file for cif_file in all_cif_files
			if cif_file.stem in allowed_pdb_ids
		]
		if not cif_files:
			raise FileNotFoundError(
				f"No CIF files found in {cif_dir} whose PDB IDs are in {pdb_id_csv}"
			)
	else:
		cif_files = all_cif_files

	# load entity -> biological map (optional). If None, will not filter by biological flag.
	# If user did not supply a path, prefer results/entity_chain_summary.csv when present.
	if entity_summary_csv is None and Path("results/entity_chain_summary.csv").exists():
		entity_summary_csv = Path("results/entity_chain_summary.csv")

	all_entity_map, biological_map = load_entity_summary(entity_summary_csv)

	file_to_accessions = {}
	file_entity_info = {}
	failed_files = []

	for cif_file in cif_files:
		try:
			extracted = extract_accessions(cif_file)
			# If we have a biological mapping for this PDB, filter to only those entity_ids
			pdb_id = cif_file.stem.upper()
			if biological_map and pdb_id in biological_map:
				allowed_entities = biological_map.get(pdb_id, set())
				filtered = [d for d in extracted if d.get("entity_id") in allowed_entities]
				# if filtering removed everything, treat as no accession
				if not filtered:
					file_to_accessions[cif_file.name] = [{"accession": "NIL", "entity_id": "NIL"}]
				else:
					file_to_accessions[cif_file.name] = filtered
			else:
				file_to_accessions[cif_file.name] = extracted

			accession_entity_ids = sorted({
				item["entity_id"]
				for item in file_to_accessions[cif_file.name]
				if item.get("entity_id") not in {None, "", "NIL"}
			})
			all_entity_ids = sorted(all_entity_map.get(pdb_id, set()))
			file_entity_info[cif_file.name] = {
				"total_entity_count": len(all_entity_ids),
				"total_entity_ids": all_entity_ids,
				"accession_entity_count": len(accession_entity_ids),
				"accession_entity_ids": accession_entity_ids,
			}
		except Exception as exc:
			failed_files.append((cif_file.name, str(exc)))
			file_to_accessions[cif_file.name] = []
			file_entity_info[cif_file.name] = {
				"total_entity_count": 0,
				"total_entity_ids": [],
				"accession_entity_count": 0,
				"accession_entity_ids": [],
			}

	# Filter out files with "NIL" accessions (no UniProt entries)
	file_to_accessions_filtered = {
		filename: accessions_dicts
		for filename, accessions_dicts in file_to_accessions.items()
		if accessions_dicts and not (len(accessions_dicts) == 1 and accessions_dicts[0].get("accession") == "NIL")
	}

	clusters = build_clusters(file_to_accessions_filtered)
	(
		json_path,
		csv_path,
		cluster_csv_path,
		homopolymer_path,
		heteropolymer_path,
		unmapped_path,
	) = write_outputs(
		output_dir,
		file_to_accessions_filtered,
		clusters,
		file_entity_info,
		file_to_accessions,
	)

	requested_count = len(cif_files)
	print(f"Parsed CIF files: {requested_count}")
	print(
		"Files with at least one UniProt accession:",
		len(file_to_accessions_filtered),
	)
	print(f"Files with NIL (no UniProt entry): {len(file_to_accessions) - len(file_to_accessions_filtered)}")
	print(f"Clusters created: {len(clusters)}")
	if allowed_pdb_ids:
		print(f"Only PDB IDs present in {pdb_id_csv} were considered for clustering.")
	print(f"Wrote: {json_path}")
	print(f"Wrote: {csv_path}")
	print(f"Wrote: {cluster_csv_path}")
	print(f"Wrote: {homopolymer_path}")
	print(f"Wrote: {heteropolymer_path}")
	print(f"Wrote: {unmapped_path}")

	if failed_files:
		print("\nFiles with parse issues (kept with empty accessions):")
		for filename, message in failed_files:
			print(f"- {filename}: {message}")


if __name__ == "__main__":
	# Collect paths interactively so the script can be used with different datasets.
	cif_dir = Path(input("Input CIF folder: ").strip())
	output_dir = Path(input("Output folder: ").strip())
	pdb_id_input = input(
		"Optional PDB ID CSV (press Enter to process all CIF files): "
	).strip()
	entity_summary_input = input(
		"Optional entity summary CSV (press Enter to skip biological filtering): "
	).strip()

	pdb_id_csv = Path(pdb_id_input) if pdb_id_input else None
	entity_summary_csv = Path(entity_summary_input) if entity_summary_input else None

	run(
		cif_dir,
		output_dir,
		pdb_id_csv=pdb_id_csv,
		entity_summary_csv=entity_summary_csv,
	)
