#!/usr/bin/env python3

import re
import csv
from pathlib import Path

import gemmi

# -------------------------------
# FASTA HEADER PARSER
# -------------------------------

def parse_fasta_header(header):
    header = header.strip()

    if not header.startswith(">"):
        return None

    header = header[1:]
    parts = header.split("|")
    if len(parts) < 3:
        return None

    # ---- PDB ID and Entity ID ----
    pdb_entity = parts[0]  # e.g. 1IWP_3
    if "_" not in pdb_entity:
        return None

    pdb_id, entity_id = pdb_entity.split("_", 1)
    if not entity_id.isdigit():
        return None

    entity_id = int(entity_id)

    # ---- Chains ----
    chains_part = parts[1].strip()
    if not (chains_part.startswith("Chain") or chains_part.startswith("Chains")):
        return None

    chains_str = (
        chains_part
        .replace("Chains", "", 1)
        .replace("Chain", "", 1)
        .strip()
    )

    label_ids = []
    auth_ids = []

    for chain_block in chains_str.split(","):
        chain_block = chain_block.strip()

        match = re.match(
            r'^([A-Za-z0-9]+)\s*(?:\[\s*auth\s+([A-Za-z0-9]+)\s*\])?$',
            chain_block
        )
        if not match:
            continue

        label = match.group(1)
        auth = match.group(2) if match.group(2) else label

        label_ids.append(label)
        auth_ids.append(auth)

    # ---- Description ----
    description = parts[2].strip()

    return {
        "pdb_id": pdb_id,
        "entity_id": entity_id,
        "label_ids": ",".join(label_ids),
        "auth_ids": ",".join(auth_ids),
        "description": description
    }


# -------------------------------
# mmCIF ENTITY POLY DETAILS
# -------------------------------

def get_entity_poly_details(block, target_entity_id):
    entity_poly = block.get_mmcif_category("_entity_poly.")

    details = {
        "is_polypeptide": False,
        "nstd_monomer": "no"
    }

    if "entity_id" not in entity_poly:
        return details

    for i, eid in enumerate(entity_poly["entity_id"]):
        if eid == str(target_entity_id):
            poly_type = entity_poly.get(
                "type", [""] * len(entity_poly["entity_id"])
            )[i]

            details["is_polypeptide"] = "polypeptide" in poly_type.lower()

            details["nstd_monomer"] = entity_poly.get(
                "nstd_monomer", ["no"] * len(entity_poly["entity_id"])
            )[i]
            break

    return details


# -------------------------------
# MAIN PROCESSING FUNCTION
# -------------------------------

def _gather_fasta_files(fasta_input):
    """Return a sorted list of FASTA files from a path (file or directory)."""

    fasta_path = Path(fasta_input)
    if fasta_path.is_file():
        return [fasta_path]

    if fasta_path.is_dir():
        fasta_files = sorted(fasta_path.glob("*.fasta"))
        if not fasta_files:
            raise FileNotFoundError(f"No FASTA files found in {fasta_path}")
        return fasta_files

    raise FileNotFoundError(f"FASTA path does not exist: {fasta_path}")


def _build_cif_lookup(cif_input):
    """Map PDB IDs (uppercased) to mmCIF paths for quick lookup."""

    cif_path = Path(cif_input)

    def _add_entry(lookup, path):
        stem = path.stem.upper()
        lookup[stem] = path
        if len(stem) > 4:
            lookup.setdefault(stem[:4], path)

    if cif_path.is_file():
        lookup = {}
        _add_entry(lookup, cif_path)
        return lookup

    if cif_path.is_dir():
        lookup = {}
        for cif_file in cif_path.rglob("*.cif"):
            _add_entry(lookup, cif_file)
        if not lookup:
            raise FileNotFoundError(f"No mmCIF files found in {cif_path}")
        return lookup

    raise FileNotFoundError(f"mmCIF path does not exist: {cif_path}")


def process_fasta_and_mmcif(fasta_input, cif_input, output_csv):
    fasta_files = _gather_fasta_files(fasta_input)
    cif_lookup = _build_cif_lookup(cif_input)
    cif_block_cache = {}

    results = []

    for fasta_file in fasta_files:
        with open(fasta_file) as f:
            for line in f:
                if not line.startswith(">"):
                    continue

                fasta_info = parse_fasta_header(line)
                if not fasta_info:
                    continue

                entity_id = fasta_info["entity_id"]
                pdb_key = fasta_info["pdb_id"].upper()

                cif_path = cif_lookup.get(pdb_key)
                if not cif_path:
                    print(f"Warning: No mmCIF file found for PDB {pdb_key}")
                    continue

                if cif_path not in cif_block_cache:
                    doc = gemmi.cif.read_file(str(cif_path))
                    cif_block_cache[cif_path] = doc.sole_block()

                block = cif_block_cache[cif_path]
                poly_info = get_entity_poly_details(block, entity_id)
                desc_lower = fasta_info["description"].lower()
                inhibitory_flag = (
                    "inhibitory protein" in desc_lower
                    or "inhibitory protien" in desc_lower
                )
                is_biological = poly_info["is_polypeptide"] and not inhibitory_flag

                record = {
                    "pdb_id": fasta_info["pdb_id"],
                    "entity_id": entity_id,
                    "label_ids": fasta_info["label_ids"],
                    "auth_ids": fasta_info["auth_ids"],
                    "description": fasta_info["description"],
                    "is_polypeptide": poly_info["is_polypeptide"],
                    "nstd_monomer": poly_info["nstd_monomer"],
                    "is_biological": is_biological
                }

                results.append(record)

    # ---- Write CSV ----
    fieldnames = [
        "pdb_id",
        "entity_id",
        "label_ids",
        "auth_ids",
        "nstd_monomer",
        "is_polypeptide",
        "is_biological",
        "description",
    ]

    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    print(f"CSV written to: {output_csv}")


# -------------------------------
# USAGE
# -------------------------------

if __name__ == "__main__": 
    fasta_file = input("Enter FASTA file path: ").strip()
    cif_file = input("Enter mmCIF file path: ").strip()
    output_csv = "entity_chain_summary.csv"

    process_fasta_and_mmcif(fasta_file, cif_file, output_csv)
