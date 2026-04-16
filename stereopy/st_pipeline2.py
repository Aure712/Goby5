#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import stereo as st
import warnings
warnings.filterwarnings('ignore')

# ====== Auto-save and unique-save utilities ======
import os
import matplotlib
matplotlib.use('Agg')  # headless friendly for matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PLOT_DIR = "plots_pdf"   # PDFs/HTMLs
TSV_DIR = "."            # keep TSVs in current dir
NPZ_DIR = "."            # keep NPZ in current dir

def _sanitize(name: str) -> str:
    return "".join(c if (c.isalnum() or c in ("-", "_", ".")) else "_" for c in str(name))

def make_unique_path(base_name: str, ext: str, outdir: str = ".") -> str:
    os.makedirs(outdir, exist_ok=True)
    base = _sanitize(base_name) if base_name else "output"
    path = os.path.join(outdir, f"{base}.{ext}")
    if not os.path.exists(path):
        return path
    i = 1
    while True:
        cand = os.path.join(outdir, f"{base}_{i}.{ext}")
        if not os.path.exists(cand):
            return cand
        i += 1

def save_current_figure(base_name: str, fig=None, ext: str = "pdf", dpi: int = 300):
    if fig is None:
        fig = plt.gcf()
    outpath = make_unique_path(base_name, ext=ext, outdir=PLOT_DIR)
    fig.savefig(outpath, dpi=dpi, bbox_inches="tight")
    print(f"[savefig] {outpath}")
    return outpath

# ========== Minimal scatter-size post-processing (方案2) ==========
def _infer_plot_kind(name: str) -> str:
    """Roughly infer plot kind by filename for selective processing."""
    low = (name or "").lower()
    if "umap" in low:
        return "umap"
    if "cluster_scatter" in low:
        return "cluster_scatter"
    if "cells_plotting" in low:
        return "cells_plotting"
    return "other"

def _auto_shrink_scatter_markers(ax, factor: float = 0.25, min_size: float = 0):
    """
    Scale down scatter marker sizes on an Axes.
    - factor: scale factor for existing sizes (e.g., 0.3 means 30%)
    - min_size: lower bound in points^2 to avoid disappearing
    """
    for coll in getattr(ax, "collections", []):
        # only PathCollection (scatter-like) has sizes
        if hasattr(coll, "get_sizes") and hasattr(coll, "set_sizes"):
            sizes = coll.get_sizes()
            if sizes is None or len(sizes) == 0:
                continue
            # sizes are in points^2
            new_sizes = np.maximum(sizes * factor, min_size)
            coll.set_sizes(new_sizes)
            # for huge scatter sets, rasterize to keep PDF smaller
            try:
                offs = coll.get_offsets() if hasattr(coll, "get_offsets") else None
                if offs is not None and len(offs) > 100000:
                    coll.set_rasterized(True)
            except Exception:
                pass

def shrink_scatter_in_current_figure(name: str, factor: float = 0.25, min_size: float = 0):
    """
    If the plot name suggests a cluster/UMAP/cells_plotting figure,
    shrink all scatter markers in the current figure.
    """
    kind = _infer_plot_kind(name)
    if kind == "other":
        return
    fig = plt.gcf()
    if fig is None:
        return
    for ax in fig.axes:
        _auto_shrink_scatter_markers(ax, factor=factor, min_size=min_size)

def step(name, func, *args, **kwargs):
    """Run a plotting function and save the current Matplotlib figure as unique PDF."""
    ret = func(*args, **kwargs)
    try:
        # shrink scatter sizes only for cluster-related plots (keep original style otherwise)
        try:
            shrink_scatter_in_current_figure(name, factor=0.3, min_size=0.5)
        except Exception as _e:
            print(f"[warn] shrink scatter failed for '{name}': {_e}")
        save_current_figure(name, ext="pdf")
    except Exception as e:
        print(f"[warn] failed to save '{name}': {e}")
    finally:
        plt.close('all')
    return ret

def save_interactive_html(obj, name: str):
    """Save interactive plot to unique HTML, best-effort for Bokeh/Panel objects."""
    try:
        from bokeh.io import output_file, save as bokeh_save
        out_html = make_unique_path(name, ext="html", outdir=PLOT_DIR)
        output_file(out_html)
        try:
            bokeh_save(obj)
        except Exception:
            layout = getattr(obj, "layout", obj)
            bokeh_save(layout)
        print(f"[savehtml] {out_html}")
        return out_html
    except Exception as e:
        print(f"[warn] could not save interactive plot '{name}' as HTML: {e}")
        return None

def save_df_tsv(df, base_name: str, sep: str = "\t", index=True):
    """Save DataFrame as TSV with unique filename in TSV_DIR."""
    path = make_unique_path(base_name, ext="tsv", outdir=TSV_DIR)
    df.to_csv(path, sep=sep, index=index)
    print(f"[savetsv] {path}")
    return path

def save_series_tsv(series, base_name: str, sep: str = "\t", index=True, header=False):
    """Save Series as TSV with unique filename in TSV_DIR."""
    path = make_unique_path(base_name, ext="tsv", outdir=TSV_DIR)
    series.to_csv(path, sep=sep, index=index, header=header)
    print(f"[savetsv] {path}")
    return path

def save_npz_unique(sparse_mat, base_name: str):
    """Save sparse matrix as NPZ with unique filename in NPZ_DIR."""
    from scipy.sparse import save_npz
    path = make_unique_path(base_name, ext="npz", outdir=NPZ_DIR)
    save_npz(path, sparse_mat)
    print(f"[savenpz] {path}")
    return path

# --- Add 'rest_mean_count' for vs.rest marker tables ---
def _parse_vs_rest_key(sKey: str, cells_df: pd.DataFrame):
    """
    Parse 'X:Y.vs.rest' style keys.
    Returns (group_key, case_value) or None if not vs.rest.
    """
    if ".vs." not in sKey:
        return None
    left, right = sKey.split(".vs.", 1)
    if right != "rest":
        return None
    if ":" in left:
        group_key, case_val = left.split(":", 1)
        return group_key, case_val
    for gk in ("DE_1", "leiden", "cluster"):
        if gk in cells_df.columns:
            return gk, left
    return None

def _compute_rest_mean_for_genes(data, group_key: str, case_val: str, genes_list, exclude_others: bool):
    """
    Compute mean raw counts over 'rest' cells for a set of genes.
    - use data.tl.raw.exp_matrix if available; else data.exp_matrix
    - exclude_others: drop rows where group_key is 'others'/'Other'
    Returns: pd.Series(index=genes_list, values=mean_count_in_rest)
    """
    cells_df = data.cells.to_df()
    if group_key not in cells_df.columns:
        return pd.Series([np.nan] * len(genes_list), index=genes_list)

    # Prefer raw counts since marker was computed with use_raw=True
    if hasattr(data.tl, "raw") and data.tl.raw is not None:
        expr_matrix = data.tl.raw.exp_matrix
        if hasattr(data.tl.raw, 'genes') and hasattr(data.tl.raw.genes, 'gene_name'):
            all_genes = data.tl.raw.genes.gene_name.tolist()
        else:
            all_genes = data.genes.gene_name.tolist()
    else:
        expr_matrix = data.exp_matrix
        all_genes = data.genes.gene_name.tolist()

    gene_to_idx = {g: i for i, g in enumerate(all_genes)}
    gene_idx = [gene_to_idx[g] for g in genes_list if g in gene_to_idx]
    valid_genes = [g for g in genes_list if g in gene_to_idx]

    col = cells_df[group_key].astype(str)
    case_val = str(case_val)
    others_vals = {"others", "Other"}
    rest_mask = (col != case_val)
    if exclude_others:
        rest_mask = rest_mask & (~col.isin(others_vals))

    rest_idx = np.where(rest_mask.values)[0]
    if rest_idx.size == 0 or len(valid_genes) == 0:
        return pd.Series([np.nan] * len(genes_list), index=genes_list)

    from scipy.sparse import issparse
    if issparse(expr_matrix):
        sub = expr_matrix.tocsr()[rest_idx, :]
        sub_g = sub[:, gene_idx]
        sums = np.array(sub_g.sum(axis=0)).ravel()
        means = sums / float(rest_idx.size)
    else:
        sub = expr_matrix[rest_idx, :]
        sub_g = sub[:, gene_idx]
        means = sub_g.mean(axis=0)

    s = pd.Series(means, index=valid_genes)
    out = pd.Series([np.nan] * len(genes_list), index=genes_list, dtype=float)
    out.loc[valid_genes] = s.values
    return out

def add_rest_mean_count_column(df_in: pd.DataFrame, sKey: str, data, exclude_others: bool):
    """
    If sKey is a '...vs.rest' key, compute and add 'rest_mean_count' column to df.
    Returns a new DataFrame.
    """
    cells_df = data.cells.to_df()
    parsed = _parse_vs_rest_key(sKey, cells_df)
    if parsed is None:
        return df_in  # not vs.rest
    group_key, case_val = parsed

    df = df_in.copy()
    if 'genes' in df.columns:
        genes_list = df['genes'].astype(str).tolist()
        rest_means = _compute_rest_mean_for_genes(data, group_key, case_val, genes_list, exclude_others)
        df['rest_mean_count'] = rest_means.values
    else:
        genes_list = df.index.astype(str).tolist()
        rest_means = _compute_rest_mean_for_genes(data, group_key, case_val, genes_list, exclude_others)
        df['rest_mean_count'] = rest_means.values

    return df
# ================================================

def main():
    print("Stereo version:", st.__version__)

    data_path = '/fast3/group_crf/home/g21shaoy23/goby5/steromics5/Tridentiger_barbel/outs/feature_expression/D03254E512.adjusted.cellbin.gef'
    st.io.read_gef_info(data_path)

    data = st.io.read_gef(file_path=data_path, bin_type='cell_bins')

    # add DE_1 from lasso
    datalasso = st.io.read_h5ad('/fast3/group_crf/home/g21shaoy23/goby5/steromics5/deg/lasso2/out/D03254E512.cellbin_1.0.h5ad', spatial_key='spatial')
    data.cells['DE_1'] = datalasso.adata.obs['DE_1']

    # QC
    data.tl.cal_qc()
    step("violin", data.plt.violin)
    step("spatial_scatter", data.plt.spatial_scatter)
    step("genes_count", data.plt.genes_count)

    # filtering
    data.tl.filter_cells(min_counts=100, min_genes=10, max_genes=10000, pct_counts_mt=5, inplace=True)

    data.tl.raw_checkpoint()  # checkpoint raw counts for downstream "use_raw"

    # normalize + log + HVGs
    data.tl.normalize_total()
    data.tl.log1p()
    data.tl.highly_variable_genes(min_mean=0.0, max_mean=5.0, min_disp=0.5, n_top_genes=10000, res_key='highly_variable_genes')
    step("highly_variable_genes", data.plt.highly_variable_genes, res_key='highly_variable_genes')

    # PCA/Neighbors/UMAP/Leiden
    data.tl.scale()
    data.tl.pca(use_highly_genes=True, n_pcs=30, res_key='pca')
    step("elbow_pca", data.plt.elbow, pca_res_key='pca')
    data.tl.neighbors(pca_res_key='pca', res_key='neighbors')
    data.tl.umap(pca_res_key='pca', neighbors_res_key='neighbors', res_key='umap')
    data.tl.leiden(neighbors_res_key='neighbors', res_key='leiden')

    step("cluster_scatter_leiden", data.plt.cluster_scatter, res_key='leiden')
    step("umap_leiden", data.plt.umap, res_key='umap', cluster_key='leiden')

    step(
        "cells_plotting_leiden",
        data.plt.cells_plotting,
        color_by='cluster',
        color_key='leiden',
        base_image='/fast3/group_crf/home/g21shaoy23/goby5/steromics5/Tridentiger_barbel/outs/image/D03254E512_ssDNA_regist.tif',
        fg_alpha=0.8
    )

    # marker: leiden vs rest
    data.tl.find_marker_genes(cluster_res_key='leiden', method='wilcoxon_test', use_highly_genes=False, use_raw=True)

    step("marker_genes_text", data.plt.marker_genes_text, res_key='marker_genes', markers_num=10, sort_key='scores')
    step("marker_genes_scatter_top5", data.plt.marker_genes_scatter, res_key='marker_genes', markers_num=5)

    _anno1 = data.plt.interact_annotation_cluster(res_cluster_key='leiden', res_marker_gene_key='marker_genes', res_key='leiden_annotation')
    save_interactive_html(_anno1, "interact_annotation_cluster")

    # filter out "others" region and redo
    celllistInclude = data.cells['DE_1'][data.cells['DE_1'] != "others"].keys()
    data.tl.filter_cells(cell_list=celllistInclude, inplace=True)

    data.tl.find_marker_genes(cluster_res_key='DE_1', method='wilcoxon_test', use_highly_genes=False, use_raw=True)
    data.tl.leiden(neighbors_res_key='neighbors', res_key='leiden')

    step("cluster_scatter_leiden_after_filter", data.plt.cluster_scatter, res_key='leiden')
    step("umap_leiden_after_filter", data.plt.umap, res_key='umap', cluster_key='leiden')
    step(
        "cells_plotting_leiden_after_filter",
        data.plt.cells_plotting,
        color_by='cluster',
        color_key='leiden',
        base_image='/fast3/group_crf/home/g21shaoy23/goby5/steromics5/Tridentiger_barbel/outs/image/D03254E512_ssDNA_regist.tif',
        fg_alpha=0.8
    )
    step("marker_genes_scatter_top5_after_filter", data.plt.marker_genes_scatter, res_key='marker_genes', markers_num=5)

    # Save DE_1 comparison marker tables (noOther) with rest_mean_count excluding 'others'
    if "marker_genes" in data.tl.result:
        for sKey in data.tl.result["marker_genes"].keys():
            if ".vs." in sKey:
                print("[marker] noOther save:", sKey)
                df_out = data.tl.result["marker_genes"][sKey]
                df_out = add_rest_mean_count_column(df_out, sKey, data, exclude_others=True)
                save_df_tsv(df_out, f"{sKey}.noOther")

    # cluster counts per DE_1 region
    cellclustersInRegions = data.cells.to_df().groupby(["DE_1", "leiden"])["leiden"].agg(["count"])
    save_df_tsv(cellclustersInRegions, "cellclustersInRegions")

    # Helper for pairwise compare (kept from original)
    def fnCompareTwoGroups(sType, sGroup1, sGroup2, data_local):
        import copy
        data2 = copy.deepcopy(data_local)
        celllistInclude2 = data2.cells[sType][(data2.cells[sType] == sGroup1) | (data2.cells[sType] == sGroup2)].keys()
        data2.tl.filter_cells(cell_list=celllistInclude2, inplace=True)

        data2.tl.find_marker_genes(
            cluster_res_key=sType, method='wilcoxon_test', use_highly_genes=False, use_raw=True,
            case_groups=sGroup1, control_groups=sGroup2
        )
        if "marker_genes" in data2.tl.result:
            for _k in data2.tl.result["marker_genes"].keys():
                if ".vs." in _k:
                    print("[marker] pairwise save:", _k)
                    df_out2 = data2.tl.result["marker_genes"][_k]
                    df_out2 = add_rest_mean_count_column(df_out2, _k, data2, exclude_others=True)
                    save_df_tsv(df_out2, f"{_k}.noOther")

    # recompute leiden markers and plots
    data.tl.leiden(neighbors_res_key='neighbors', res_key='leiden')
    step("cluster_scatter_leiden_again", data.plt.cluster_scatter, res_key='leiden')
    step("umap_leiden_again", data.plt.umap, res_key='umap', cluster_key='leiden')

    # groups 1..20
    for g in map(str, range(1, 21)):
        step(f"cluster_scatter_leiden_group_{g}", data.plt.cluster_scatter, res_key='leiden', groups=[g])

    step(
        "cells_plotting_leiden_again",
        data.plt.cells_plotting,
        color_by='cluster',
        color_key='leiden',
        base_image='/fast3/group_crf/home/g21shaoy23/goby5/steromics5/Tridentiger_barbel/outs/image/D03254E512_ssDNA_regist.tif',
        fg_alpha=0.8
    )

    data.tl.find_marker_genes(cluster_res_key='leiden', method='wilcoxon_test', use_highly_genes=False, use_raw=True)
    step("marker_genes_text_again", data.plt.marker_genes_text, res_key='marker_genes', markers_num=10, sort_key='scores')
    step("marker_genes_scatter_top5_again", data.plt.marker_genes_scatter, res_key='marker_genes', markers_num=5)

    _anno2 = data.plt.interact_annotation_cluster(res_cluster_key='leiden', res_marker_gene_key='marker_genes', res_key='leiden_annotation')
    save_interactive_html(_anno2, "interact_annotation_cluster_again")

    ins = data.plt.interact_spatial_scatter(width=1000, height=800, poly_select=True)
    save_interactive_html(ins, "interact_spatial_scatter")
    # do not call ins.show() in batch mode

    # Save latest marker results (likely for leiden) with rest_mean_count (not excluding 'others' here)
    if "marker_genes" in data.tl.result:
        for sKey in data.tl.result["marker_genes"].keys():
            if ".vs." in sKey:
                print("[marker] save:", sKey)
                df_out = data.tl.result["marker_genes"][sKey]
                df_out = add_rest_mean_count_column(df_out, sKey, data, exclude_others=False)
                save_df_tsv(df_out, f"{sKey}")

    # ============ Marker genes total expression extraction ============
    def extract_marker_genes_total_expression(data_local, output_prefix="marker_genes"):
        all_marker_genes = set()

        if "marker_genes" not in data_local.tl.result:
            print("No marker_genes result found in data.tl.result")
            return None

        marker_keys = data_local.tl.result["marker_genes"].keys()
        print(f"Found {len(marker_keys)} marker gene comparisons")

        for key in marker_keys:
            if ".vs." in key:
                marker_df = data_local.tl.result["marker_genes"][key]
                if 'genes' in marker_df.columns:
                    genes = marker_df['genes'].tolist()
                else:
                    genes = marker_df.index.tolist()
                all_marker_genes.update(genes)
                print(f"{key}: {len(genes)} genes")

        all_marker_genes = sorted(list(all_marker_genes))
        print(f"\nTotal unique marker genes: {len(all_marker_genes)}")

        if hasattr(data_local.tl, 'raw') and data_local.tl.raw is not None:
            print("Raw data checkpoint exists, using raw data for full gene set")
            expr_matrix = data_local.tl.raw.exp_matrix
            if hasattr(data_local.tl.raw, 'genes'):
                gene_names = data_local.tl.raw.genes.gene_name.tolist()
            elif hasattr(data_local.tl.raw, 'gene_names'):
                gene_names = data_local.tl.raw.gene_names.tolist()
            else:
                gene_names = data_local.genes.gene_name.tolist()
        else:
            print("Using current expression matrix")
            expr_matrix = data_local.exp_matrix
            gene_names = data_local.genes.gene_name.tolist()

        n_cells = expr_matrix.shape[0]
        n_genes = expr_matrix.shape[1]
        print(f"\nExpression matrix shape: (cells={n_cells}, genes={n_genes})")
        print(f"Number of gene names available: {len(gene_names)}")

        if len(gene_names) != n_genes:
            print(f"WARNING: Gene names count ({len(gene_names)}) doesn't match matrix columns ({n_genes})")
            if len(gene_names) > n_genes:
                print("This might happen after highly_variable_genes filtering")
                print("Attempting to match gene names to matrix...")
            else:
                print("Error: Fewer gene names than matrix columns")
                return None

        gene_to_idx = {gene: idx for idx, gene in enumerate(gene_names)}
        gene_indices = []
        genes_found = []
        genes_not_found = []

        for gene in all_marker_genes:
            if gene in gene_to_idx:
                idx = gene_to_idx[gene]
                if idx < n_genes:
                    gene_indices.append(idx)
                    genes_found.append(gene)
                else:
                    genes_not_found.append(gene)
            else:
                genes_not_found.append(gene)

        print(f"\nGenes found in expression matrix: {len(genes_found)}")
        if genes_not_found:
            print(f"Genes not found: {len(genes_not_found)}")
            if len(genes_not_found) <= 20:
                print(f"Not found: {genes_not_found}")
            else:
                print(f"First 20 not found: {genes_not_found[:20]}")

        if genes_found:
            print(f"\nCalculating total expression for {len(genes_found)} marker genes...")

            if hasattr(expr_matrix, 'tocsc'):
                expr_matrix_csc = expr_matrix.tocsc()
                marker_expr = expr_matrix_csc[:, gene_indices]
                total_expression = np.array(marker_expr.sum(axis=0)).flatten()
                mean_expression = np.array(marker_expr.mean(axis=0)).flatten()
                n_cells_expressed = np.array((marker_expr > 0).sum(axis=0)).flatten()
                if hasattr(marker_expr, 'toarray'):
                    marker_expr_dense = marker_expr.toarray()
                    max_expression = marker_expr_dense.max(axis=0)
                else:
                    max_expression = np.array(marker_expr.max(axis=0)).flatten()
            else:
                marker_expr = expr_matrix[:, gene_indices]
                total_expression = marker_expr.sum(axis=0)
                mean_expression = marker_expr.mean(axis=0)
                max_expression = marker_expr.max(axis=0)
                n_cells_expressed = (marker_expr > 0).sum(axis=0)

            if len(total_expression.shape) > 1:
                total_expression = total_expression.flatten()
            if len(mean_expression.shape) > 1:
                mean_expression = mean_expression.flatten()
            if len(max_expression.shape) > 1:
                max_expression = max_expression.flatten()
            if len(n_cells_expressed.shape) > 1:
                n_cells_expressed = n_cells_expressed.flatten()

            pct_cells_expressed = (n_cells_expressed / n_cells) * 100

            expression_summary = pd.DataFrame({
                'gene': genes_found,
                'total_expression': total_expression,
                'mean_expression': mean_expression,
                'max_expression': max_expression,
                'n_cells_expressed': n_cells_expressed,
                'pct_cells_expressed': pct_cells_expressed,
                'n_total_cells': n_cells
            }).sort_values('total_expression', ascending=False)

            save_df_tsv(expression_summary, f"{output_prefix}_total")

            print(f"\nTop 10 genes by total expression:")
            top10 = expression_summary.head(10)[['gene', 'total_expression', 'mean_expression', 'pct_cells_expressed']]
            for _, row in top10.iterrows():
                print(f"  {row['gene']}: total={row['total_expression']:.2f}, mean={row['mean_expression']:.4f}, pct_cells={row['pct_cells_expressed']:.2f}%")

            try:
                cells_df = data_local.cells.to_df()

                if 'leiden' in cells_df.columns:
                    print("\nCalculating total expression by leiden clusters...")
                    cluster_expression = {}
                    if hasattr(expr_matrix, 'tocsc'):
                        expr_matrix_csc = expr_matrix.tocsc()
                    for cluster in sorted(cells_df['leiden'].unique()):
                        cluster_cells = cells_df[cells_df['leiden'] == cluster].index
                        cell_indices = [i for i, cell_id in enumerate(cells_df.index) if cell_id in cluster_cells]
                        if len(cell_indices) > 0:
                            if hasattr(expr_matrix, 'tocsc'):
                                cluster_expr = expr_matrix_csc[cell_indices, :][:, gene_indices]
                                cluster_total = np.array(cluster_expr.sum(axis=0)).flatten()
                                cluster_mean = np.array(cluster_expr.mean(axis=0)).flatten()
                            else:
                                cluster_expr = expr_matrix[cell_indices, :][:, gene_indices]
                                cluster_total = cluster_expr.sum(axis=0)
                                cluster_mean = cluster_expr.mean(axis=0)
                            cluster_expression[f"cluster_{cluster}_total"] = cluster_total
                            cluster_expression[f"cluster_{cluster}_mean"] = cluster_mean
                            cluster_expression[f"cluster_{cluster}_n_cells"] = len(cell_indices)

                    if cluster_expression:
                        cluster_df = pd.DataFrame(cluster_expression, index=genes_found)
                        save_df_tsv(cluster_df, f"{output_prefix}_by_cluster")

                if 'DE_1' in cells_df.columns:
                    print("\nCalculating total expression by DE_1 regions...")
                    region_expression = {}
                    region_expression['all_regions_total'] = total_expression
                    region_expression['all_regions_mean'] = mean_expression
                    region_expression['all_regions_n_cells'] = n_cells

                    if hasattr(expr_matrix, 'tocsc'):
                        expr_matrix_csc = expr_matrix.tocsc()

                    regions = [r for r in cells_df['DE_1'].unique() if r != 'others']
                    for region in regions:
                        region_cells = cells_df[cells_df['DE_1'] == region].index
                        cell_indices = [i for i, cell_id in enumerate(cells_df.index) if cell_id in region_cells]
                        if len(cell_indices) > 0:
                            if hasattr(expr_matrix, 'tocsc'):
                                region_expr = expr_matrix_csc[cell_indices, :][:, gene_indices]
                                region_total = np.array(region_expr.sum(axis=0)).flatten()
                                region_mean = np.array(region_expr.mean(axis=0)).flatten()
                            else:
                                region_expr = expr_matrix[cell_indices, :][:, gene_indices]
                                region_total = region_expr.sum(axis=0)
                                region_mean = region_expr.mean(axis=0)
                            region_expression[f"{region}_total"] = region_total
                            region_expression[f"{region}_mean"] = region_mean
                            region_expression[f"{region}_n_cells"] = len(cell_indices)

                    if 'others' in cells_df['DE_1'].unique():
                        others_cells = cells_df[cells_df['DE_1'] == 'others'].index
                        cell_indices = [i for i, cell_id in enumerate(cells_df.index) if cell_id in others_cells]
                        if len(cell_indices) > 0:
                            if hasattr(expr_matrix, 'tocsc'):
                                others_expr = expr_matrix_csc[cell_indices, :][:, gene_indices]
                                others_total = np.array(others_expr.sum(axis=0)).flatten()
                                others_mean = np.array(others_expr.mean(axis=0)).flatten()
                            else:
                                others_expr = expr_matrix[cell_indices, :][:, gene_indices]
                                others_total = others_expr.sum(axis=0)
                                others_mean = others_expr.mean(axis=0)
                            region_expression['others_total'] = others_total
                            region_expression['others_mean'] = others_mean
                            region_expression['others_n_cells'] = len(cell_indices)

                    if region_expression:
                        region_df = pd.DataFrame(region_expression, index=genes_found)
                        cols = region_df.columns.tolist()
                        all_regions_cols = [col for col in cols if col.startswith('all_regions')]
                        other_cols = [col for col in cols if not col.startswith('all_regions')]
                        region_df = region_df[all_regions_cols + other_cols]
                        save_df_tsv(region_df, f"{output_prefix}_by_region")

                        print(f"  Regions included: {', '.join(regions)}")
                        if 'others' in cells_df['DE_1'].unique():
                            print(f"  Also included 'others' region")
                        print(f"  Total expression column: 'all_regions_total'")

            except Exception as e:
                print(f"Could not calculate cluster/region expression: {e}")

            return expression_summary

        else:
            print("No marker genes found in expression matrix")
            return None

    def verify_data_structure(data_local):
        print("\n" + "="*60)
        print("DATA STRUCTURE VERIFICATION")
        print("="*60)
        print(f"\nExpression matrix shape: {data_local.exp_matrix.shape}")

        cells_df = data_local.cells.to_df()
        print(f"Cells DataFrame shape: {cells_df.shape}")
        print(f"Number of cells: {cells_df.shape[0]}")
        print(f"Number of gene names: {len(data_local.genes.gene_name)}")

        if data_local.exp_matrix.shape[0] == cells_df.shape[0]:
            print("\n✓ Confirmed: Matrix rows = cells")
        if data_local.exp_matrix.shape[1] == len(data_local.genes.gene_name):
            print("✓ Confirmed: Matrix columns = genes")

        print(f"\nMatrix format: ({data_local.exp_matrix.shape[0]} cells × {data_local.exp_matrix.shape[1]} genes)")

        if hasattr(data_local.tl, 'raw') and data_local.tl.raw is not None:
            print(f"\nRaw data available:")
            print(f"  Raw matrix shape: {data_local.tl.raw.exp_matrix.shape}")
            if hasattr(data_local.tl.raw, 'genes'):
                print(f"  Raw gene count: {len(data_local.tl.raw.genes.gene_name)}")

        if 'highly_variable_genes' in data_local.tl.result:
            print("\n⚠ Note: highly_variable_genes filtering has been applied")
            hvg_result = data_local.tl.result['highly_variable_genes']
            if hasattr(hvg_result, 'highly_variable'):
                try:
                    n_hvg = hvg_result.highly_variable.sum()
                    print(f"  Number of highly variable genes: {n_hvg}")
                except Exception:
                    pass
        print("="*60 + "\n")

    def extract_all_genes_total_expression(data_local, output_prefix="all_genes"):
        print("\nExtracting total expression for ALL genes...")

        expr_matrix = data_local.exp_matrix
        gene_names = data_local.genes.gene_name.tolist()

        n_cells = expr_matrix.shape[0]
        n_genes = expr_matrix.shape[1]

        print(f"Processing {n_genes} genes across {n_cells} cells...")

        if hasattr(expr_matrix, 'tocsc'):
            expr_matrix_csc = expr_matrix.tocsc()
            total_expr = np.array(expr_matrix_csc.sum(axis=0)).flatten()
            mean_expr = np.array(expr_matrix_csc.mean(axis=0)).flatten()
            n_cells_expr = np.array((expr_matrix_csc > 0).sum(axis=0)).flatten()
        else:
            total_expr = expr_matrix.sum(axis=0)
            mean_expr = expr_matrix.mean(axis=0)
            n_cells_expr = (expr_matrix > 0).sum(axis=0)

        if len(total_expr.shape) > 1:
            total_expr = total_expr.flatten()
        if len(mean_expr.shape) > 1:
            mean_expr = mean_expr.flatten()
        if len(n_cells_expr.shape) > 1:
            n_cells_expr = n_cells_expr.flatten()

        all_genes_df = pd.DataFrame({
            'gene': gene_names[:n_genes],
            'total_expression': total_expr,
            'mean_expression': mean_expr,
            'n_cells_expressed': n_cells_expr,
            'pct_cells_expressed': (n_cells_expr / n_cells) * 100
        }).sort_values('total_expression', ascending=False)

        save_df_tsv(all_genes_df, f"{output_prefix}_total_expression")
        print(f"Total genes processed: {len(all_genes_df)}")

        return all_genes_df

    # Main execution
    print("=" * 60)
    print("MARKER GENES TOTAL EXPRESSION EXTRACTION")
    print("=" * 60)

    verify_data_structure(data)
    marker_expr_summary = extract_marker_genes_total_expression(data, output_prefix="marker_genes")

    if marker_expr_summary is not None:
        print("\n" + "=" * 60)
        print("SUCCESS: Marker genes extraction complete!")
        print(f"Extracted total expression for {len(marker_expr_summary)} marker genes")
        print("Files created (unique names):")
        print("  - marker_genes_total*.tsv (main results)")
        print("  - marker_genes_by_cluster*.tsv (if leiden clustering exists)")
        print("  - marker_genes_by_region*.tsv (if DE_1 regions exist)")
        print("=" * 60)
    else:
        print("\n" + "=" * 60)
        print("WARNING: Some marker genes could not be found")
        print("Extracting ALL genes as alternative...")
        print("=" * 60)
        all_genes_expr = extract_all_genes_total_expression(data, 'all_genes')
        if all_genes_expr is not None:
            print(f"\nExtracted expression for {len(all_genes_expr)} total genes")
            print("You can filter 'all_genes_total_expression*.tsv' for your marker genes")

    # Per-cell and per-cluster counts exports
    from scipy.sparse import issparse, csr_matrix

    X = data.exp_matrix
    genes = list(data.genes.gene_name)
    cells = data.cells.to_df().index.astype(str).tolist()

    gene_total = np.asarray(X.sum(axis=0)).ravel() if issparse(X) else X.sum(axis=0)
    save_df_tsv(pd.DataFrame({'gene': genes, 'total_expression': gene_total}), "current_region_gene_total", index=False)

    X_to_save = X.tocsr() if issparse(X) else csr_matrix(X)
    save_npz_unique(X_to_save, "per_cell_all_expression")
    save_series_tsv(pd.Series(cells), "per_cell_ids", index=False, header=False)
    save_series_tsv(pd.Series(genes), "gene_names", index=False, header=False)
    print('Saved: current_region_gene_total*.tsv, per_cell_all_expression*.npz, per_cell_ids*.tsv, gene_names*.tsv')

    # raw counts per region/cluster
    raw_mat = data.tl.raw.exp_matrix
    gene_name = data.tl.raw.genes.gene_name.tolist()

    if issparse(raw_mat):
        gene_total_raw = np.ravel(raw_mat.sum(axis=0))
    else:
        gene_total_raw = raw_mat.sum(axis=0)

    save_df_tsv(pd.DataFrame({'gene': gene_name, 'total_counts_in_region': gene_total_raw}), "region_total_counts", index=False)
    print('① 已保存（唯一命名）：region_total_counts*.tsv')

    cells_df = data.cells.to_df()
    leiden = cells_df['leiden'].astype(str)
    cluster_gene_counts = pd.DataFrame(index=gene_name, columns=sorted(leiden.unique()), dtype=np.int64)

    for cl in cluster_gene_counts.columns:
        mask = leiden == cl
        if mask.any():
            sub_mat = raw_mat[mask.values, :]
            if issparse(sub_mat):
                cluster_gene_counts[cl] = np.ravel(sub_mat.sum(axis=0))
            else:
                cluster_gene_counts[cl] = sub_mat.sum(axis=0)

    n_cells_by_cluster = cells_df.groupby('leiden').size().reindex(cluster_gene_counts.columns).fillna(0)
    cluster_gene_counts.loc['__n_cells__'] = n_cells_by_cluster.values

    save_df_tsv(cluster_gene_counts, "cluster_gene_counts")
    print('② 已保存（唯一命名）：cluster_gene_counts*.tsv  (行=基因, 列=leiden 聚类, 最后一行=__n_cells__)')

if __name__ == "__main__":
    main()
