#!/usr/bin/env python3
import sys
import os
import csv
import gemmi
from typing import Optional, Dict, Any, List

def extract_chain_sequence(cif_path: str, chain_id: str, entity_id: str) -> Optional[Dict[str, Any]]:
    # Basic input checks
    if not chain_id or not chain_id.strip():
        print("Error: chain_id cannot be empty")
        return None
    if not entity_id or not entity_id.strip():
        print("Error: entity_id cannot be empty")
        return None

    target_chain  = str(chain_id).strip()
    target_entity = str(entity_id).strip()

    # Read CIF
    try:
        doc = gemmi.cif.read_file(cif_path)
        block = doc.sole_block()
    except Exception as e:
        print(f"Error reading mmCIF {cif_path}: {e}")
        return None

    # Find _atom_site loop via reliable tags
    col = block.find_loop('_atom_site.label_entity_id')
    if not col:
        col = block.find_loop('_atom_site.auth_asym_id')
        if not col:
            print("Error: _atom_site loop not found via expected tags")
            return None

    loop = col.get_loop()
    if loop is None or loop.length() == 0:
        print("Error: _atom_site loop is missing or empty")
        return None

    # Build name→index map from loop.tags
    tags: List[str] = list(loop.tags)
    tag_to_idx = {t: i for i, t in enumerate(tags)}

    def first_index(candidates: List[str]) -> int:
        for t in candidates:
            if t in tag_to_idx:
                return tag_to_idx[t]
        return -1

    # Resolve indices with robust fallbacks
    idx_entity_id = first_index(['_atom_site.label_entity_id', '_atom_site.auth_entity_id'])
    idx_asym_id   = first_index(['_atom_site.label_asym_id',   '_atom_site.auth_asym_id'])
    idx_seq_id    = first_index(['_atom_site.auth_seq_id',     '_atom_site.label_seq_id'])
    idx_group_pdb = first_index(['_atom_site.group_PDB'])
    idx_comp_id   = first_index(['_atom_site.label_comp_id',   '_atom_site.auth_comp_id'])
    idx_ins_code  = first_index(['_atom_site.label_ins_code',  '_atom_site.pdbx_PDB_ins_code'])

    # Validate required columns
    required = [idx_entity_id, idx_asym_id, idx_seq_id, idx_group_pdb, idx_comp_id]
    if any(i < 0 for i in required):
        print("Error: Missing one or more required _atom_site columns.")
        return None

    # Iterate rows with numeric indexing and tuple access
    found = False
    observed_positions = set()
    modified_residue_map: Dict[int, str] = {}
    has_insertion_code = False

    nrows = loop.length()
    ncols = len(loop.tags)
    values = list(loop.values)
    for r in range(nrows):
        base = r * ncols
        if base + ncols > len(values):
            break
        if values[base + idx_entity_id] == target_entity and values[base + idx_asym_id] == target_chain:
            # seq_id may be '.' for solvent or non-polymer; guard int conversion
            try:
                pos = int(values[base + idx_seq_id])
            except ValueError:
                continue
            group_pdb = values[base + idx_group_pdb]
            resname   = values[base + idx_comp_id]
            ins_code  = values[base + idx_ins_code] if idx_ins_code >= 0 else '.'

            observed_positions.add(pos)
            if group_pdb == 'HETATM' and resname != 'HOH':
                modified_residue_map[pos] = resname
            if ins_code and ins_code not in ['.', '?']:
                has_insertion_code = True

            found = True

    if not found:
        print(f"No atoms found for chain '{target_chain}' and entity '{target_entity}'")
        return None

    # Summaries
    obs_sorted = sorted(observed_positions)
    modified_residues = [{"three_letter": n, "position": p}
                         for p, n in sorted(modified_residue_map.items())]

    no_internal_missing = False
    if obs_sorted:
        expected = list(range(obs_sorted[0], obs_sorted[-1] + 1))
        no_internal_missing = (obs_sorted == expected)

    return {
        "observed_positions": obs_sorted,
        "non_standard_residues": modified_residues,
        "no_internal_missing": no_internal_missing,
        "has_insertion_code": has_insertion_code
    }

def parse_bool(value: str) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}

def parse_label_ids(value: str) -> List[str]:
    if value is None:
        return []
    raw = str(value).strip().strip('"')
    if not raw:
        return []
    return [part.strip() for part in raw.split(",") if part.strip()]

def find_cif_path(cif_dir: str, pdb_id: str) -> Optional[str]:
    if not pdb_id:
        return None
    candidates = [
        f"{pdb_id}.cif",
        f"{pdb_id.lower()}.cif",
        f"{pdb_id.upper()}.cif",
        f"{pdb_id}.cif.gz",
        f"{pdb_id.lower()}.cif.gz",
        f"{pdb_id.upper()}.cif.gz",
    ]
    for name in candidates:
        path = os.path.join(cif_dir, name)
        if os.path.exists(path):
            return path
    return None

def main():
    if len(sys.argv) < 3:
        print("Usage: script.py <path/to/summary.csv> <path/to/cif_dir> [output.csv]")
        print("Example: script.py results/entity_chain_summary.csv data/unzipped_cif results/gemi_output.csv")
        csv_path = input("Enter path to summary CSV: ").strip()
        cif_dir = input("Enter path to CIF directory: ").strip()
        out_path = input("Enter output CSV path (leave blank for gemi_output.csv): ").strip()
        if not out_path:
            out_path = "gemi_output.csv"
    else:
        csv_path = sys.argv[1]
        cif_dir = sys.argv[2]
        out_path = sys.argv[3] if len(sys.argv) >= 4 else "gemi_output.csv"
    
    # Create log file path
    if out_path.endswith(os.sep) or os.path.isdir(out_path):
        log_path = os.path.join(out_path, "gemi_processing.log")
        out_path = os.path.join(out_path, "gemi_output.csv")
    else:
        out_dir = os.path.dirname(out_path) if os.path.dirname(out_path) else "."
        out_basename = os.path.splitext(os.path.basename(out_path))[0]
        log_path = os.path.join(out_dir, f"{out_basename}.log")

    if not os.path.exists(csv_path):
        print(f"Error: CSV file not found: {csv_path}")
        sys.exit(1)
    if not os.path.isdir(cif_dir):
        print(f"Error: CIF directory not found: {cif_dir}")
        sys.exit(1)

    # Create log directory if needed
    log_dir = os.path.dirname(log_path)
    if log_dir and not os.path.exists(log_dir):
        os.makedirs(log_dir, exist_ok=True)
    
    results: List[Dict[str, Any]] = []
    log_messages: List[str] = []
    with open(csv_path, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if not parse_bool(row.get("is_biological")):
                continue

            pdb_id = str(row.get("pdb_id", "")).strip()
            entity_id = str(row.get("entity_id", "")).strip()
            chains = parse_label_ids(row.get("label_ids"))

            if not pdb_id or not entity_id or not chains:
                continue

            cif_path = find_cif_path(cif_dir, pdb_id)
            if not cif_path:
                log_messages.append(f"WARNING: CIF not found for {pdb_id} in {cif_dir}")
                for chain_id in chains:
                    log_messages.append(f"  Skipped: {pdb_id} entity={entity_id} chain={chain_id} (cif_not_found)")
                continue

            log_messages.append(f"Processing {pdb_id}: Using CIF file {cif_path}")
            
            for chain_id in chains:
                out = extract_chain_sequence(cif_path, chain_id, entity_id)
                if out is None:
                    log_messages.append(f"  WARNING: No atoms found for {pdb_id} entity={entity_id} chain={chain_id}")
                    continue

                obs_pos = out["observed_positions"]
                start_position = str(obs_pos[0]) if obs_pos else ""
                end_position = str(obs_pos[-1]) if obs_pos else ""
                
                non_standard = ";".join(
                    f"{r['three_letter']}:{r['position']}" for r in out["non_standard_residues"]
                )
                results.append({
                    "pdb_id": pdb_id,
                    "entity_id": entity_id,
                    "chain_id": chain_id,
                    "start_position": start_position,
                    "end_position": end_position,
                    "non_standard_residues": non_standard,
                    "no_internal_missing": str(out["no_internal_missing"]),
                    "has_insertion_code": str(out["has_insertion_code"]),
                })
                log_messages.append(f"  Successfully processed: {pdb_id} entity={entity_id} chain={chain_id} (positions {start_position}-{end_position})")

    # Write log file
    with open(log_path, "w") as log_file:
        log_file.write("\n".join(log_messages))
    
    # Write CSV output
    out_dir = os.path.dirname(out_path)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    fieldnames = [
        "pdb_id",
        "entity_id",
        "chain_id",
        "start_position",
        "end_position",
        "non_standard_residues",
        "no_internal_missing",
        "has_insertion_code",
    ]

    with open(out_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    print(f"Wrote {len(results)} records to {out_path}")
    print(f"Wrote log to {log_path}")

if __name__ == "__main__":
    main()
