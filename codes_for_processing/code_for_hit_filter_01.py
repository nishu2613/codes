import os
import csv
import subprocess
from Bio import SeqIO
from Bio.Blast import NCBIXML

class BLASTProcessor:
    def __init__(self, xml_folder, blast_db, output_folder, query_fasta, no_hits_log):
        self.xml_folder = xml_folder        
        self.blast_db = blast_db            
        self.output_folder = output_folder  
        self.query_fasta = query_fasta      
        self.no_hits_log = no_hits_log      
        self.csv_output_folder = os.path.join(self.output_folder, "filtered_hits_csv")

        os.makedirs(self.output_folder, exist_ok=True)
        os.makedirs(self.csv_output_folder, exist_ok=True)
        
        no_hits_log_dir = os.path.dirname(self.no_hits_log)
        if no_hits_log_dir:
            os.makedirs(no_hits_log_dir, exist_ok=True)
        
        with open(self.no_hits_log, "w") as log:
            log.write("XML_File\tQuery_Definition\n")
        
        self.no_hits_count = 0  

    def get_hit_metrics(self, hit, query_length):
        """
        Evaluates a hit holistically across all its HSPs to solve alignment fragmentation.
        Calculates true non-overlapping coverage, weighted identity, and minimum E-value.
        """
        # 1. Merge overlapping query intervals for exact non-redundant coverage
        intervals = sorted(
            (min(h.query_start, h.query_end), max(h.query_start, h.query_end))
            for h in hit.hsps
        )
        merged_intervals = [intervals[0]]
        for start, end in intervals[1:]:
            prev_start, prev_end = merged_intervals[-1]
            if start <= prev_end + 1:
                merged_intervals[-1] = (prev_start, max(prev_end, end))
            else:
                merged_intervals.append((start, end))

        covered_query_bases = sum(end - start + 1 for start, end in merged_intervals)
        cumulative_coverage = covered_query_bases / query_length

        # 2. Calculate overall weighted identity and lowest E-value
        total_identities = sum(h.identities for h in hit.hsps)
        total_align_length = sum(h.align_length for h in hit.hsps)
        weighted_identity = total_identities / total_align_length if total_align_length > 0 else 0
        
        # The E-value acts as the guardrail against fragmented noise
        best_evalue = min(h.expect for h in hit.hsps)

        return weighted_identity, cumulative_coverage, best_evalue

    def process_all_xml_files(self):
        xml_files = [f for f in os.listdir(self.xml_folder) if f.endswith(".xml")]
        if not xml_files:
            print("No XML files found.")
            return

        for i, xml_file in enumerate(xml_files, 1):
            print(f"\n => Processing file {i}/{len(xml_files)}: {xml_file}")
            self.process_xml_file(xml_file)

        print("\n" + "="*60)
        print(f" All BLAST XML files processed!")
        print(f" Total processed: {len(xml_files)}")
        print(f" Queries with hits: {len(xml_files) - self.no_hits_count}")
        print(f" Queries with NO hits: {self.no_hits_count}")
        print(f" Filtered outputs saved to: {self.output_folder}")
        print(f" No-hits log saved to: {self.no_hits_log}")
        print("="*60)

    def process_xml_file(self, xml_file):
        xml_path = os.path.join(self.xml_folder, xml_file)
        with open(xml_path, "r") as file:
            blast_record = NCBIXML.read(file)

        query_def = blast_record.query
        query_length = blast_record.query_length
        hit_details = []

        for hit in blast_record.alignments:
            if not hit.hsps:
                continue

            # Apply hit-level aggregation logic
            identity, coverage, evalue = self.get_hit_metrics(hit, query_length)

            # Strict pan-genomic filtration thresholds
            if identity >= 0.90 and coverage >= 0.70 and evalue < 1e-5:
                clean_id = hit.hit_id.strip()
                hit_details.append((clean_id, identity, coverage, evalue))

        hit_def_list = sorted(set(hit[0] for hit in hit_details))
        
        if not hit_def_list:
            print(f"   └─ ⚠ No hits found for this query. Skipping file creation.")
            with open(self.no_hits_log, "a") as log:
                log.write(f"{xml_file}\t{query_def}\n")
            self.no_hits_count += 1
            return
        
        print(f"   Found {len(hit_def_list)} valid strain-specific ortholog hits")
        self.extract_sequences_fast(query_def, hit_def_list, xml_file)
        self.write_hit_defs_to_csv(hit_details, xml_file)

    def extract_sequences_fast(self, query_def, hit_def_list, xml_file):
        """Uses blastdbcmd to extract sequences and cleanly restores original headers on the fly."""
        fasta_id = os.path.splitext(xml_file)[0]
        query_fasta_file = os.path.join(self.query_fasta, f"{fasta_id}.fasta")
        
        if not os.path.exists(query_fasta_file):
            print(f"   └─ ❌ ERROR: Query FASTA file not found! Checked path: {query_fasta_file}")
            return
        
        query_id = None
        query_sequence = None
        with open(query_fasta_file, "r") as query_file:
            for record in SeqIO.parse(query_file, "fasta"):
                if query_def and query_def in record.description:
                    query_id = record.id
                    query_sequence = str(record.seq)
                    break

        if not query_id:
            return

        safe_query_id = query_id.replace("|", "_")
        output_file = os.path.join(self.output_folder, f"filtered_hits_{safe_query_id}.faa")

        # 1. Create temporary batch file for blastdbcmd
        temp_batch_file = os.path.join(self.output_folder, f"temp_{safe_query_id}_batch.txt")
        with open(temp_batch_file, "w") as batch_f:
            for hit_id in hit_def_list:
                batch_f.write(f"{hit_id}\n")

        # 2. Write the Query Sequence first
        with open(output_file, "w") as outfile:
            outfile.write(f">{query_def}\n")
            for i in range(0, len(query_sequence), 60):
                outfile.write(query_sequence[i:i+60] + "\n")

        # 3. Run blastdbcmd, capture text in Python, clean up headers, and write out
        cmd = [
            "blastdbcmd",
            "-db", self.blast_db,
            "-entry_batch", temp_batch_file,
            "-outfmt", "%f"
        ]
        
        print(f"   ├─ 🚀 Querying BLAST index and purifying headers...")
        
        # Run process and capture text stdout stream and errors
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.stderr:
            print(f"   └─ ⚠ blastdbcmd warning/error: {result.stderr.strip()}")
        
        with open(output_file, "a") as outfile:
            for line in result.stdout.splitlines():
                if line.startswith(">"):
                    outfile.write(f"{line}\n")
                else:
                    outfile.write(f"{line}\n")

        # 4. Clean up the temporary batch file
        if os.path.exists(temp_batch_file):
            os.remove(temp_batch_file)

    def write_hit_defs_to_csv(self, hit_details, xml_file):
        if not hit_details:
            return
        base_name = os.path.splitext(xml_file)[0]
        csv_file = os.path.join(self.csv_output_folder, f"{base_name}_hit_defs.csv")
        with open(csv_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["Hit_Def", "Identity", "Coverage", "E-value"])
            for hit_def, identity, coverage, evalue in sorted(hit_details):
                writer.writerow([hit_def, f"{identity:.4f}", f"{coverage:.4f}", f"{evalue:.2e}"])

if __name__ == "__main__":
    xml_folder = input("Enter the path to the folder containing BLAST XML files: ").strip()
    blast_db = input("Enter the path to the reformatted BLAST database: ").strip()
    output_folder = input("Enter the path where you want to save the output files: ").strip()
    query_fasta = input("Enter the path to the folder containing query FASTA files: ").strip()
    no_hits_log = input("Enter the path for the no-hits log file: ").strip()

    processor = BLASTProcessor(xml_folder, blast_db, output_folder, query_fasta, no_hits_log)
    processor.process_all_xml_files()