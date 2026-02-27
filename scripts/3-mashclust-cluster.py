#!/usr/bin/env python3
"""
Script 3: Cluster Genomes and Select Representatives (Multi-Reference Support)
"""
import sys
import argparse
import json
import numpy as np
import re
from pathlib import Path
from collections import defaultdict

def extract_genome_id(path_str):
    """Extract clean ID (GCF/GCA) or fallback to folder name"""
    path = Path(path_str)
    # Regex changed: [._] matches either dot or underscore before version
    match = re.search(r'(GC[FA]_\d{9}[._]\d+)', path.name)
    if match: return match.group(1)
    # Attempt 2: Parent folder name
    match = re.search(r'(GC[FA]_\d{9}[._]\d+)', path.parent.name)
    if match: return match.group(1)
    # Fallback
    return path.parent.name

def normalize_id(acc_id):
    """Normalize Accession ID by replacing dots with underscores."""
    return acc_id.replace('.', '_')

def main():
    parser = argparse.ArgumentParser(description="Step 3: Cluster Genomes")
    parser.add_argument("input_dir", help="Directory with distances.txt or distances.txt file")
    parser.add_argument("-o", "--output", required=True, help="Output file (representatives.txt)")
    parser.add_argument("-t", "--threshold", type=float, default=0.9997)
    parser.add_argument("-n", "--num-representatives", type=int, default=5)
    parser.add_argument("--reference", nargs='+', help="Reference Accession(s) to protect")
    parser.add_argument("--no-reference-protection", action="store_true")
    args = parser.parse_args()

    # If -o is /path/to/representatives.txt, we write to that file and put json in that folder
    out_file = Path(args.output)
    out_dir = out_file.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Identify distances.txt location
    input_path = Path(args.input_dir)
    dist_file = input_path if input_path.is_file() else input_path / "distances.txt"

    if not dist_file.exists():
        print(f"[ERROR] Distance file not found at: {dist_file}")
        sys.exit(1)

    print(f"[INFO] Threshold: {args.threshold} | Max Reps: {args.num_representatives}")
    
    # 1. Parse Distance Matrix
    neighbors = defaultdict(set)
    genome_paths = {} # Map ID -> Full Path
    
    distance_threshold = 1.0 - args.threshold
    
    print("[INFO] Parsing distance matrix...")
    with open(dist_file, 'r') as f:
        # Read header
        header = f.readline().strip().split('\t')[1:]
        all_genomes_list = [p.strip() for p in header]
        
        # Map IDs to paths
        for p in all_genomes_list:
            gid = extract_genome_id(p)
            genome_paths[gid] = p
            
        genome_ids_ordered = list(genome_paths.keys())

        # Read rows
        for i, line in enumerate(f):
            parts = line.strip().split('\t')
            if len(parts) < 2: continue
            
            query_id = extract_genome_id(parts[0])
            dists = parts[1:]
            
            for j, d_str in enumerate(dists):
                try:
                    d = float(d_str)
                    target_id = genome_ids_ordered[j]
                    
                    if query_id != target_id and d <= distance_threshold:
                        neighbors[query_id].add(target_id)
                        neighbors[target_id].add(query_id)
                except ValueError:
                    continue

    # 2. Identify References (Improved with normalization)
    ref_ids = set()
    if not args.no_reference_protection and args.reference:
        for ref_input in args.reference:
            ref_input_clean = normalize_id(ref_input)
            found = False
            for gid in genome_paths.keys():
                if ref_input_clean in normalize_id(gid):
                    ref_ids.add(gid)
                    print(f"[INFO] Reference genome identified: {gid} (matched {ref_input})")
                    found = True
                    break 
            if not found:
                print(f"[WARNING] Reference {ref_input} not found in dataset.")

    # 3. Greedy Clustering
    print("[INFO] Running Greedy Clustering...")
    clusters = []
    assigned = set()
    
    # Sort genomes: References first (forced), then by connectivity
    genome_connectivity = [(g, len(neighbors[g])) for g in genome_paths.keys()]
    
    if ref_ids:
        refs = [x for x in genome_connectivity if x[0] in ref_ids]
        others = [x for x in genome_connectivity if x[0] not in ref_ids]
        refs.sort(key=lambda x: x[1], reverse=True)
        others.sort(key=lambda x: x[1], reverse=True)
        genome_connectivity = refs + others
        print(f"[INFO] {len(refs)} references prioritized in clustering queue.")
    else:
        genome_connectivity.sort(key=lambda x: x[1], reverse=True)

    for g, _ in genome_connectivity:
        if g in assigned:
            continue
            
        cluster = {g}
        assigned.add(g)
        
        for neighbor in neighbors[g]:
            if neighbor not in assigned:
                cluster.add(neighbor)
                assigned.add(neighbor)
        
        clusters.append(cluster)

    print(f"[INFO] Created {len(clusters)} clusters.")

    # 4. Select Representatives
    representatives = []
    np.random.seed(42) # Reproducibility
    
    for cluster in clusters:
        cluster_list = list(cluster)
        
        if len(cluster_list) <= args.num_representatives:
            representatives.extend(cluster_list)
        else:
            selected = []
            pool = list(cluster_list)
            
            if ref_ids:
                refs_in_cluster = [m for m in pool if m in ref_ids]
                for r in refs_in_cluster:
                    selected.append(r)
                    pool.remove(r)
                    if len(selected) >= args.num_representatives:
                        break
            
            slots_left = args.num_representatives - len(selected)
            if slots_left > 0:
                picked = np.random.choice(pool, size=slots_left, replace=False)
                selected.extend(picked)
            
            representatives.extend(selected)

    # 5. Output
    # A) Representatives list (Writes to the exact file Snakemake wants)
    with open(out_file, "w") as f:
        for r in representatives:
            f.write(f"{genome_paths[r]}\n")
            
    # B) JSON data (In the same folder as representatives.txt)
    json_path = out_dir / "clustering_data.json"
    json_data = {
        'clusters': [list(c) for c in clusters],
        'neighbors': {k: list(v) for k, v in neighbors.items()},
        'genome_names': genome_paths,
        'identity_threshold': args.threshold
    }
    
    with open(json_path, "w") as f:
        json.dump(json_data, f)
        
    print(f"[SUCCESS] Selected {len(representatives)} representatives.")

if __name__ == "__main__":
    main()