#!/usr/bin/env python3
"""
Script 1: Sketch and Filter (Target vs Non-Target)
"""
import subprocess
import sys
import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Step 1: Sketch and Filter")
    parser.add_argument("genomes_dir", help="Directory containing genome folders from Step 0")
    parser.add_argument("-o", "--output-dir", required=True)
    parser.add_argument("-f", "--filter", required=True, help="Bacteria name prefix (e.g., mycobacterium_tuberculosis)")
    parser.add_argument("-k", "--kmer", type=int, default=31)
    parser.add_argument("-s", "--sketch-size", type=int, default=100000)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    target_genomes = []
    non_target_genomes = []
    
    print(f"[INFO] Scanning genomes in {args.genomes_dir}")
    print(f"[INFO] Filtering for folders starting with: '{args.filter}'")

    base_path = Path(args.genomes_dir)
    
    for genome_file in base_path.rglob("GENOME_*.fna"):
        # Parent folder name (e.g., mycobacterium_tuberculosis_1773_GCF...)
        parent_folder = genome_file.parent.name.lower()
        
        if parent_folder.startswith(args.filter.lower()):
            target_genomes.append(genome_file)
        else:
            non_target_genomes.append(genome_file)

    print(f"[INFO] Targets found (to be sketched): {len(target_genomes)}")
    print(f"[INFO] Non-targets found (to be moved directly to results): {len(non_target_genomes)}")

    # Save lists for later steps
    with open(out_dir / "targets.txt", "w") as f:
        for p in target_genomes:
            f.write(f"{p}\n")
            
    with open(out_dir / "non_targets.txt", "w") as f:
        for p in non_target_genomes:
            f.write(f"{p}\n")

    if not target_genomes:
        print("[ERROR] No target genomes found! Check your --filter parameter or accessions.")
        sys.exit(1)

    # Mash Sketch ONLY on Targets
    sketch_file = out_dir / "genomes_sketch"
    
    cmd = [
        "mash", "sketch",
        "-s", str(args.sketch_size),
        "-k", str(args.kmer),
        "-p", str(args.threads),
        "-o", str(sketch_file),
        "-l", str(out_dir / "targets.txt")
    ]

    print("[INFO] Running Mash Sketch on targets...")
    try:
        subprocess.run(cmd, check=True)
        print(f"[SUCCESS] Sketch created: {sketch_file}.msh")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Mash sketch failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()