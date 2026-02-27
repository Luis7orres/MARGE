#!/usr/bin/env python3
"""
Script 2: Calculate Distances (Zero-RAM Streaming Version)
"""
import subprocess
import sys
import argparse
import time
from pathlib import Path

def main():
    # 1. Argument Parsing
    parser = argparse.ArgumentParser(description="Step 2: Calculate Distances")
    parser.add_argument("input_dir", help="Directory with .msh file")
    parser.add_argument("-o", "--output-dir", required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--chunk-size", type=int, default=2000, help="Compatibility flag")
    args = parser.parse_args()

    # 2. Setup Paths
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    input_dir = Path(args.input_dir)

    sketch_path = input_dir / "genomes_sketch.msh"
    dist_file = out_dir / "distances.txt"

    if not sketch_path.exists():
        print(f"[ERROR] Sketch file not found: {sketch_path}")
        sys.exit(1)

    # 3. Define Mash Command
    # We pass the sketch twice to trigger All-vs-All comparison
    cmd = [
        "mash", "dist",
        "-t",                   # Tabular output (square matrix)
        "-p", str(args.threads),
        str(sketch_path),       # Reference set
        str(sketch_path)        # Query set
    ]

    print("-" * 60)
    print(f"[INFO] Starting calculation (All-vs-All)")
    print(f"[INFO] Threads: {args.threads}")
    print(f"[INFO] Output:  {dist_file}")
    print(f"[INFO] Logic:   Streaming directly to disk to bypass Python RAM limits.")
    print("-" * 60)

    # 4. Execution via Direct Streaming
    start_time = time.time()
    try:
        # We open the file and pass the handle directly to the subprocess.
        # This ensures Mash writes to the disk and Python never sees the data.
        with open(dist_file, "w") as f_out:
            result = subprocess.run(
                cmd, 
                stdout=f_out, 
                stderr=subprocess.PIPE, 
                text=True,
                check=True
            )
            
        duration = time.time() - start_time
        print(f"[SUCCESS] Calculation finished in {duration/60:.2f} minutes.")
        print(f"[SUCCESS] Matrix saved to {dist_file}")

    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Mash dist failed.")
        print(f"[ERROR] Stderr output: {e.stderr}")
        sys.exit(1)
    except Exception as e:
        print(f"[ERROR] Unexpected error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()