#!/usr/bin/env python3
"""
Script 4: Visualización de Heatmap Clusterizado con Escala de Color Específica
Límite Azul: 0.06 | Transición gradual a 0.07 (Amarillo) | Rojo: 0.1
Mantiene TODAS las funcionalidades originales (Quicktree, MDS, Plots).
"""
import json
import argparse
import sys
import subprocess
from pathlib import Path
import numpy as np

# Intentar importaciones de visualización
try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
    print("[WARNING] Plotly no encontrado.")

# Intentar importaciones de análisis
try:
    from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram
    from scipy.spatial.distance import squareform
    from sklearn.manifold import MDS
    ANALYSIS_LIBS_AVAILABLE = True
except ImportError:
    ANALYSIS_LIBS_AVAILABLE = False
    print("[WARNING] Scipy o Scikit-learn no encontrados.")

def extract_genome_name(path_str):
    return Path(path_str).parent.name

def generate_quicktree(output_dir, dist_file, representatives_file):
    if not Path(representatives_file).exists(): return
    rep_names = [extract_genome_name(line.strip()) for line in open(representatives_file) if line.strip()]
    rep_set = set(rep_names)

    dist_path = Path(dist_file)
    genome_order = []
    distance_matrix = []

    with open(dist_path, 'r') as f:
        header = f.readline().strip().split('\t')[1:]
        rep_indices = [i for i, path in enumerate(header) if extract_genome_name(path) in rep_set]
        for line in f:
            parts = line.strip().split('\t')
            if extract_genome_name(parts[0]) in rep_set:
                genome_order.append(extract_genome_name(parts[0]))
                distance_matrix.append([float(parts[i+1]) for i in rep_indices])

    phylip_file = output_dir / "representatives_phylip.dist"
    n = len(genome_order)
    with open(phylip_file, 'w') as f:
        f.write(f"{n}\n")
        for i, name in enumerate(genome_order):
            dists_str = '\t'.join(f"{d:.6f}" for d in distance_matrix[i])
            f.write(f"{name}\t{dists_str}\n")

    tree_file = output_dir / "representatives_tree.nwk"
    try:
        cmd = ["quicktree", "-in", "m", "-out", "t", str(phylip_file)]
        with open(tree_file, 'w') as f_out:
            subprocess.run(cmd, check=True, stdout=f_out, stderr=subprocess.PIPE)
    except Exception: pass

def create_plots(data, output_dir):
    if not PLOTLY_AVAILABLE: return
    clusters = data['clusters']
    cluster_sizes = [len(c) for c in clusters]
    fig = make_subplots(rows=1, cols=2, subplot_titles=('Tamaño de Clústeres', 'Acumulado'))
    fig.add_trace(go.Histogram(x=cluster_sizes, marker_color='steelblue'), row=1, col=1)
    sorted_sizes = sorted(cluster_sizes, reverse=True)
    fig.add_trace(go.Scatter(y=np.cumsum(sorted_sizes), mode='lines'), row=1, col=2)
    fig.write_html(output_dir / "cluster_distribution.html")

def generate_advanced_visualizations(output_dir, dist_file, representatives_file):
    if not PLOTLY_AVAILABLE or not ANALYSIS_LIBS_AVAILABLE: return

    # 1. Cargar y filtrar matriz de distancias
    rep_set = set(extract_genome_name(line.strip()) for line in open(representatives_file) if line.strip())
    with open(dist_file, 'r') as f:
        header = f.readline().strip().split('\t')[1:]
        indices = [i for i, h in enumerate(header) if extract_genome_name(h) in rep_set]
        labels = [extract_genome_name(header[i]) for i in indices]
        data = []
        for line in f:
            parts = line.strip().split('\t')
            if extract_genome_name(parts[0]) in rep_set:
                data.append([float(parts[i+1]) for i in indices])
    
    matrix_np = np.array(data)
    n = len(labels)

    # 2. Clustering Jerárquico (UPGMA)
    link = linkage(squareform(matrix_np, checks=False), method='average')
    order = leaves_list(link)
    ordered_labels = [labels[i] for i in order]
    matrix_ord = matrix_np[order, :][:, order]

    # --- ESCALA DE COLORES MEJORADA (transición más visible, menos blanco) ---
    cs = [
        [0.00, '#313695'], # 0.00: Azul oscuro
        [0.40, '#74add1'], # 0.04: Azul medio
        [0.60, '#abd9e9'], # 0.06: Azul muy claro (94% ANI)
        [0.65, '#e0f3f8'], # 0.065: Azul muy pálido (NO blanco puro)
        [0.70, '#fee08b'], # 0.07: Amarillo claro (transición visible)
        [0.80, '#fdae61'], # 0.08: Naranja
        [0.90, '#f46d43'], # 0.09: Naranja-Rojo
        [1.00, '#d73027']  # 0.10: Rojo
    ]

    # --- 3) Preparar dendrogramas SciPy (siguiendo el método que funciona) ---
    # Top dendrogram
    dtop = dendrogram(link, labels=ordered_labels, orientation='top', no_plot=True)
    # Right dendrogram
    dright = dendrogram(link, labels=ordered_labels, orientation='right', no_plot=True)

    # Escalas de posición de hoja de SciPy: hojas en 5,15,25,...,10*(n-1)+5
    # Para pasar a 0..n-1: pos_scaled = (pos_leaf - 5) / 10
    icoord_top = dtop['icoord']
    dcoord_top = dtop['dcoord']
    max_d = max(max(dc) for dc in dcoord_top)

    icoord_r = dright['icoord']
    dcoord_r = dright['dcoord']

    # --- 4) Crear figura ---
    fig_hm = go.Figure()
    
    # Dendrograma Superior (arriba del heatmap)
    top_traces = []
    for icoord, dcoord in zip(icoord_top, dcoord_top):
        xs = [(x - 5)/10 for x in icoord]  # escala 0..n-1
        ys = dcoord                         # altura de cluster
        top_traces.append(go.Scatter(
            x=xs, y=ys,
            mode='lines',
            line=dict(color='black', width=1.2),
            hoverinfo='none',
            showlegend=False
        ))

    # Dendrograma Derecho (a la derecha del heatmap) - CORREGIDO
    right_traces = []
    for icoord, dcoord in zip(icoord_r, dcoord_r):
        # La clave es invertir correctamente las coordenadas Y
        ys = [(n-1) - (y - 5)/10 for y in icoord]  # filas invertidas
        xs = dcoord                                 # distancia
        right_traces.append(go.Scatter(
            x=xs, y=ys,
            mode='lines',
            line=dict(color='black', width=1.2),
            hoverinfo='none',
            showlegend=False
        ))

    # Añadir dendrograma superior
    for tr in top_traces:
        fig_hm.add_trace(tr.update(xaxis='x2', yaxis='y2'))
    
    # Añadir dendrograma derecho
    for tr in right_traces:
        fig_hm.add_trace(tr.update(xaxis='x3', yaxis='y3'))

    # Heatmap central
    fig_hm.add_trace(go.Heatmap(
        z=matrix_ord, 
        x=ordered_labels, 
        y=ordered_labels,
        colorscale=cs,
        zmin=0, zmax=0.1, 
        xaxis='x', yaxis='y',
        colorbar=dict(
            title=dict(
                text="<b>Mash Distance</b>", 
                font=dict(size=17)
            ),
            tickvals=[0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1],
            ticktext=['0', '0.01', '0.02', '0.03', '0.04', '0.05', '0.06', '0.07', '0.08', '0.09', '0.1'],
            tickfont=dict(size=16),
            thickness=14,
            x=0.92, 
            len=0.84
        ),
        hovertemplate='<b>Q:</b> %{y}<br><b>R:</b> %{x}<br><b>Dist:</b> %{z:.4f}<extra></extra>'
    ))

    # --- 5) Layout con ejes correctamente alineados ---
    fig_hm.update_layout(
        title=dict(
            text="<b>Dual Clustered Mash Heatmap</b>",
            x=0.5, 
            font=dict(size=18), 
            y=0.99
        ),
        width=1400, 
        height=1200, 
        showlegend=False,
        
        # Heatmap central
        xaxis=dict(
            domain=[0.08, 0.82], 
            zeroline=False, 
            showgrid=False,
            showticklabels=False,
            range=[-0.5, n-0.5]
        ),
        yaxis=dict(
            domain=[0.08, 0.92], 
            zeroline=False, 
            showgrid=False,
            autorange='reversed',
            showticklabels=False,
            range=[-0.5, n-0.5]
        ),
        
        # Dendrograma Superior
        xaxis2=dict(
            domain=[0.08, 0.82], 
            anchor='y2',
            zeroline=False, 
            showticklabels=False, 
            showgrid=False,
            range=[-0.5, n-0.5]
        ),
        yaxis2=dict(
            domain=[0.92, 0.99], 
            anchor='x2',
            zeroline=False, 
            showticklabels=False, 
            showgrid=False,
            range=[0, max_d]
        ),
        
        # Dendrograma Derecho - CORREGIDO
        xaxis3=dict(
            domain=[0.82, 0.90], 
            anchor='y3',
            zeroline=False, 
            showticklabels=False, 
            showgrid=False,
            range=[0, max_d]
        ),
        yaxis3=dict(
            domain=[0.08, 0.92],  # Mismo rango que el heatmap
            anchor='x3',
            zeroline=False, 
            showticklabels=False, 
            showgrid=False,
            range=[-0.5, n-0.5]   # Mismo rango que el heatmap
        ),
        
        paper_bgcolor='white', 
        plot_bgcolor='white',
        margin=dict(l=30, r=120, t=20, b=40)
    )
    fig_hm.write_html(output_dir / "representatives_clustered_heatmap.html")

    # PCoA (MDS) - Se mantiene exactamente igual
    try:
        mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42)
        coords = mds.fit_transform(matrix_np)
        fig_pcoa = go.Figure(go.Scatter(
            x=coords[:, 0], y=coords[:, 1], mode='markers+text',
            text=labels, textposition="top center",
            marker=dict(size=12, color=coords[:, 0], colorscale='Viridis', line=dict(width=1, color='white'))
        ))
        fig_pcoa.update_layout(title="<b>PCoA of Representatives</b>", template="plotly_white", width=1100, height=800)
        fig_pcoa.write_html(output_dir / "representatives_pcoa.html")
    except Exception: pass

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", help="Directorio con clustering_data.json")
    parser.add_argument("-o", "--output-dir", required=True)
    args = parser.parse_args()

    in_dir = Path(args.input_dir)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    json_path = in_dir / "clustering_data.json"
    rep_txt = in_dir / "representatives.txt" 
    dist_txt = in_dir.parent / "2-distances" / "distances.txt"

    if not json_path.exists():
        print(f"[ERROR] No se encontró {json_path}")
        sys.exit(1)

    with open(json_path, 'r') as f:
        data = json.load(f)
        
    create_plots(data, out_dir)
    generate_quicktree(out_dir, dist_txt, rep_txt)
    generate_advanced_visualizations(out_dir, dist_txt, rep_txt)
    print(f"\n[DONE] Archivos generados en: {out_dir}")

if __name__ == "__main__":
    main()