#!/usr/bin/env python
"""
Script 5: Finalize Results (The "Perfect Integration" approach)
Usa directamente dataset-manager.py para garantizar compatibilidad 100% con Pipeline 1.
"""
import argparse
import shutil
import re
import os
import sys
import subprocess
import time
from pathlib import Path

def normalize_accession(path_str):
    """Extrae GCA_000000000.1 de cualquier cadena."""
    pattern = r'(GC[FA])_(\d{9})[_.](\d+)'
    match = re.search(pattern, path_str)
    if match:
        return f"{match.group(1)}_{match.group(2)}.{match.group(3)}"
    return None

def main():
    parser = argparse.ArgumentParser(description="Step 5: Finalize Results")
    parser.add_argument("--non-targets", required=True)
    parser.add_argument("--representatives", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--acc-file", required=True, help="Ruta para selected_accessions.txt")
    parser.add_argument("--genomes-subdir", default="uncompressed")
    # NUEVO ARGUMENTO: Ruta a tu script dataset-manager.py
    parser.add_argument("--dataset-manager", required=True, help="Ruta absoluta a dataset-manager.py")
    args = parser.parse_args()

    print("\n" + "="*60)
    print(" INICIANDO INTEGRACIÓN CON DATASET-MANAGER")
    print("="*60)

    # 1. Verificar herramientas
    datasets_exe = shutil.which("datasets")
    if not datasets_exe:
        print("[ERROR CRÍTICO] No se encontró el binario 'datasets'.")
        return

    dm_script = Path(args.dataset_manager)
    if not dm_script.exists():
        print(f"[ERROR CRÍTICO] No se encuentra dataset-manager.py en: {dm_script}")
        return

    # 2. Configurar rutas
    base_dir = Path(args.out_dir) 
    uncompressed_root = base_dir / args.genomes_subdir
    base_dir.mkdir(parents=True, exist_ok=True)
    uncompressed_root.mkdir(parents=True, exist_ok=True)

    # 3. Leer listas
    selected_items = []
    def collect(file_path):
        if not Path(file_path).exists(): return
        with open(file_path, 'r') as f:
            for line in f:
                path_str = line.strip()
                if not path_str: continue
                acc = normalize_accession(path_str)
                # folder_name es el ID de carpeta que usa tu pipeline (staphylococcus_..._GCA_...)
                folder_name = Path(path_str).parent.name if "/" in path_str else Path(path_str).name
                if acc:
                    selected_items.append((acc, folder_name))

    collect(args.non_targets)
    collect(args.representatives)

    if not selected_items:
        print("[ERROR] No hay genomas para procesar.")
        return

    print(f"[INFO] Procesando {len(selected_items)} genomas...")

    final_accessions = []
    index_lines = [] # Para el índice de kSNP
    failed = 0

    for i, (acc, folder_name) in enumerate(selected_items, 1):
        zip_dest = base_dir / f"{folder_name}.zip"
        # dataset-manager crea una carpeta basada en el nombre del zip o argumentos,
        # así que forzaremos que la salida sea uncompressed/
        
        # A. Descarga
        if not (zip_dest.exists() and zip_dest.stat().st_size > 5000):
            print(f"[{i}] Descargando {acc}...")
            # IMPORTANTE: Incluimos cds, gff3 y seq-report para que dataset-manager no falle al buscarlos
            cmd_dl = [datasets_exe, "download", "genome", "accession", acc, 
                      "--include", "genome,seq-report,cds,gff3", 
                      "--filename", str(zip_dest)]
            try:
                subprocess.run(cmd_dl, check=True, capture_output=False)
            except subprocess.CalledProcessError:
                print(f"    [!] Fallo descarga NCBI.")
                failed += 1
                continue
        else:
             print(f"[{i}] {acc} ya descargado.")

        # B. Ejecutar dataset-manager.py
        # El comando equivalente en ksnp-preparation.sh es:
        # python3 dataset-manager.py build-dataset --output <dir> <zip>
        print(f"    -> Ejecutando dataset-manager...")
        
        cmd_dm = [
            sys.executable, str(dm_script), 
            "build-dataset", 
            "--output", str(uncompressed_root),
            str(zip_dest)
        ]

        try:
            # dataset-manager imprime a stdout la línea del reporte (Genus Species Strain...)
            # Necesitamos capturarla para construir el índice
            result = subprocess.run(cmd_dm, check=True, capture_output=True, text=True)
            
            # La salida de dataset-manager es una línea tabulada.
            # ksnp-preparation.sh le añade el filename al principio.
            # Salida DM:  Genus \t Species \t Strain \t TaxID \t Path
            # Output Final: FolderName \t Genus \t Species ...
            
            dm_output = result.stdout.strip()
            
            if dm_output:
                # Construir la línea del índice kSNP
                full_index_line = f"{folder_name}\t{dm_output}"
                index_lines.append(full_index_line)
                final_accessions.append(acc)
                print("    [OK] Procesado correctamente.")
            else:
                print("    [!] dataset-manager no devolvió datos.")
                failed += 1

        except subprocess.CalledProcessError as e:
            print(f"    [!] Error en dataset-manager: {e.stderr}")
            failed += 1
        except Exception as e:
            print(f"    [!] Error inesperado: {e}")
            failed += 1

    # --- GUARDAR ARCHIVOS ---

    # 1. selected_accessions.txt
    print(f"\n[INFO] Guardando {args.acc_file}")
    with open(args.acc_file, 'w') as f:
        for acc in sorted(list(set(final_accessions))):
            f.write(f"{acc}\n")

    # 2. ksnp_files_index.tsv (Fundamental para que el Pipeline 1 funcione)
    # Lo guardamos en el mismo directorio que selected_accessions.txt
    ksnp_index_path = Path(args.acc_file).parent / "ksnp_files_index.tsv"
    print(f"[INFO] Guardando índice kSNP en {ksnp_index_path}")
    with open(ksnp_index_path, 'w') as f:
        for line in index_lines:
            f.write(f"{line}\n")

    print(f"\n[FIN] Completados: {len(final_accessions)} | Fallidos: {failed}")

if __name__ == "__main__":
    main()