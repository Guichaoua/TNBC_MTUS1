from __future__ import annotations

import json
import math
import re
from io import StringIO
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from scipy.stats import linregress, pearsonr


PROJECT = Path(__file__).resolve().parents[1]
MANUSCRIPT = PROJECT / "manuscript"
OUT = PROJECT / "figures"
OUT.mkdir(parents=True, exist_ok=True)

NOTEBOOK_DIR = PROJECT / "notebooks"
CELL_LINES = PROJECT / "data" / "raw" / "Cell-lines"


mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 8.5,
        "axes.titlesize": 10,
        "axes.labelsize": 8.5,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "legend.fontsize": 7.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.7,
        "savefig.bbox": "tight",
    }
)


CORE_TERMS = [
    "HALLMARK_MYC_TARGETS_V1",
    "HALLMARK_MYC_TARGETS_V2",
    "HALLMARK_DNA_REPAIR",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
]

TERM_LABELS = {
    "HALLMARK_MYC_TARGETS_V1": "MYC targets V1",
    "HALLMARK_MYC_TARGETS_V2": "MYC targets V2",
    "HALLMARK_DNA_REPAIR": "DNA repair",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION": "Oxidative phosphorylation",
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE": "Unfolded protein response",
    "HALLMARK_G2M_CHECKPOINT": "G2M checkpoint",
    "HALLMARK_E2F_TARGETS": "E2F targets",
}

COHORT_LABELS = {
    "SRP042620": "SRP042620",
    "VUMC": "VUMC",
    "GSE192341": "GSE192341",
    "GSE181466": "GSE181466",
    "TCGA": "TCGA",
    "GSE202203": "GSE202203",
    "SRP157974": "SRP157974",
    "TNBC_All_Pooled": "TNBC pooled",
}


def clean_columns(df: pd.DataFrame) -> pd.DataFrame:
    if isinstance(df.columns, pd.MultiIndex):
        df.columns = [c[0] for c in df.columns]
    if "Term" in df.columns:
        df = df.rename(columns={"Term": "cohort"})
    return df


def notebook_html_tables(notebook_name: str, cell_index: int) -> list[pd.DataFrame]:
    nb = json.loads((NOTEBOOK_DIR / notebook_name).read_text())
    html = ""
    for output in nb["cells"][cell_index].get("outputs", []):
        data = output.get("data", {})
        if "text/html" in data:
            value = data["text/html"]
            html += "".join(value) if isinstance(value, list) else str(value)
    if not html:
        return []
    return [clean_columns(t) for t in pd.read_html(StringIO(html))]


def extract_tnbc_gsea() -> tuple[pd.DataFrame, pd.DataFrame]:
    nes = notebook_html_tables("terciles_log2_DE_GSEA_v3.ipynb", 24)[0]
    fdr = notebook_html_tables("terciles_log2_DE_GSEA_v3.ipynb", 28)[0]
    nes = nes.set_index("cohort")
    fdr = fdr.set_index("cohort")
    nes = nes.loc[list(COHORT_LABELS.keys()), CORE_TERMS]
    fdr = fdr.loc[list(COHORT_LABELS.keys()), CORE_TERMS]
    nes.to_csv(OUT / "source_tnbc_core_gsea_nes.csv")
    fdr.to_csv(OUT / "source_tnbc_core_gsea_fdr.csv")
    return nes, fdr


def save_figure(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUT / f"{stem}.pdf")
    fig.savefig(OUT / f"{stem}.png", dpi=600)
    plt.close(fig)


def fdr_to_size(fdr: float) -> float:
    fdr = max(float(fdr), 1e-5)
    score = min(-math.log10(fdr), 5)
    return 20 + score * 38


def plot_tnbc_core_dotplot(nes: pd.DataFrame, fdr: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 3.8))
    xs, ys, colors, sizes = [], [], [], []
    for y, term in enumerate(CORE_TERMS[::-1]):
        for x, cohort in enumerate(nes.index):
            xs.append(x)
            ys.append(y)
            colors.append(nes.loc[cohort, term])
            sizes.append(fdr_to_size(fdr.loc[cohort, term]))

    sc = ax.scatter(
        xs,
        ys,
        c=colors,
        s=sizes,
        cmap="RdBu_r",
        norm=TwoSlopeNorm(vmin=-3.8, vcenter=0, vmax=3.8),
        edgecolor="#333333",
        linewidth=0.45,
    )
    ax.set_xticks(range(len(nes.index)))
    ax.set_xticklabels([COHORT_LABELS[x] for x in nes.index], rotation=45, ha="right")
    ax.set_yticks(range(len(CORE_TERMS)))
    ax.set_yticklabels([TERM_LABELS[t] for t in CORE_TERMS[::-1]])
    ax.set_title("Reproducible pathway activity in MTUS1-low TNBC tumours", loc="left", weight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.grid(axis="x", color="#e8e8e8", linewidth=0.7)
    ax.grid(axis="y", color="#efefef", linewidth=0.7)
    ax.set_axisbelow(True)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    cbar = fig.colorbar(sc, ax=ax, pad=0.015, fraction=0.035)
    cbar.set_label("NES")
    legend_fdrs = [0.05, 0.01, 0.001, 0.0001]
    handles = [
        plt.scatter([], [], s=fdr_to_size(v), facecolor="white", edgecolor="#333333", linewidth=0.45)
        for v in legend_fdrs
    ]
    labels = ["FDR <= 0.05", "FDR <= 0.01", "FDR <= 0.001", "FDR <= 0.0001"]
    ax.legend(handles, labels, title="Adjusted p-value", frameon=False, loc="upper left", bbox_to_anchor=(1.13, 1.02))
    save_figure(fig, "Figure1_TNBC_core_GSEA")


def read_depmap_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = pd.read_csv(CELL_LINES / "res" / "meta_24cellLines_chronos.csv")
    ids = meta["DepMap_ID"].dropna().astype(str).tolist()
    expr = pd.read_csv(
        CELL_LINES / "OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected.csv",
        index_col=0,
        usecols=lambda c: c in ["", "MTUS1 (57509)"] or c.startswith("Unnamed"),
    )
    if "MTUS1 (57509)" not in expr.columns:
        expr = pd.read_csv(CELL_LINES / "OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected.csv", index_col=0)
        expr = expr[["MTUS1 (57509)"]]
    effects = pd.read_csv(CELL_LINES / "CRISPRGeneEffect.csv", index_col=0)
    common = [i for i in ids if i in expr.index and i in effects.index]
    expr = expr.loc[common].rename(columns={"MTUS1 (57509)": "MTUS1_expression"})
    effects = effects.loc[common]
    return expr, effects


def compute_depmap_correlations(expr: pd.DataFrame, effects: pd.DataFrame) -> pd.DataFrame:
    x = expr["MTUS1_expression"]
    rows = []
    for col in effects.columns:
        y = effects[col]
        ok = x.notna() & y.notna()
        if ok.sum() < 10:
            continue
        r, p = pearsonr(x[ok], y[ok])
        rows.append(
            {
                "Gene": col.split(" (", 1)[0],
                "column": col,
                "Pearson Correlation": r,
                "p_value": p,
                "n": int(ok.sum()),
                "min_score": float(y[ok].min()),
                "max_score": float(y[ok].max()),
            }
        )
    corr = pd.DataFrame(rows)
    corr = corr[corr["min_score"] < -1].copy()
    corr["minus_log10_p"] = -np.log10(corr["p_value"].clip(lower=1e-300))
    selected = pd.read_csv(PROJECT / "data" / "derived" / "depmap_mtus1_associated_essential_genes.csv", sep=";")
    selected_genes = set(selected["Gene"].astype(str))
    corr["selected_candidate"] = corr["Gene"].isin(selected_genes)
    corr.to_csv(OUT / "source_depmap_all_essentiality_correlations.csv", index=False)
    selected.to_csv(OUT / "source_depmap_selected_candidates.csv", index=False)
    return corr


def plot_depmap_prioritisation(expr: pd.DataFrame, effects: pd.DataFrame, corr: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(7.4, 6.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.05], width_ratios=[1.0, 1.05], hspace=0.42, wspace=0.36)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])

    selected = corr[corr["selected_candidate"]]
    ax1.hist(corr["Pearson Correlation"], bins=42, color="#d7dee8", edgecolor="white", linewidth=0.5)
    ax1.hist(selected["Pearson Correlation"], bins=np.linspace(corr["Pearson Correlation"].min(), corr["Pearson Correlation"].max(), 42), color="#d76545", alpha=0.85)
    ax1.axvline(0, color="#777777", linewidth=0.8)
    ax1.axvline(0.4, color="#9b2f2f", linestyle="--", linewidth=0.9)
    ax1.set_title("A. Correlation screen", loc="left", weight="bold")
    ax1.set_xlabel("Pearson r: Chronos score vs MTUS1 expression")
    ax1.set_ylabel("Genes")
    ax1.text(
        0.98,
        0.93,
        f"{len(selected)} candidates",
        transform=ax1.transAxes,
        ha="right",
        va="top",
        color="#9b2f2f",
    )

    ax2.scatter(
        corr["Pearson Correlation"],
        corr["minus_log10_p"],
        s=12,
        color="#b9c3cf",
        edgecolor="none",
        alpha=0.75,
    )
    ax2.scatter(
        selected["Pearson Correlation"],
        selected["minus_log10_p"],
        s=20,
        color="#d76545",
        edgecolor="white",
        linewidth=0.3,
        alpha=0.95,
    )
    ax2.axhline(-math.log10(0.05), color="#555555", linestyle="--", linewidth=0.8)
    ax2.axvline(0, color="#777777", linewidth=0.7)
    for gene in ["MYC", "DMAP1", "DDX10", "IMP3", "MCM6", "DKC1"]:
        row = corr[corr["Gene"] == gene]
        if not row.empty:
            r = row.iloc[0]
            ax2.text(r["Pearson Correlation"] + 0.015, r["minus_log10_p"] + 0.05, gene, fontsize=7)
    ax2.set_title("B. Candidate dependency genes", loc="left", weight="bold")
    ax2.set_xlabel("Pearson r")
    ax2.set_ylabel("-log10(p)")

    myc_col = "MYC (4609)"
    dep = -effects.loc[expr.index, myc_col]
    x = expr["MTUS1_expression"]
    ok = x.notna() & dep.notna()
    slope, intercept, r_value, p_value, _, = linregress(x[ok], dep[ok])
    xx = np.linspace(float(x[ok].min()), float(x[ok].max()), 100)
    ax3.scatter(x[ok], dep[ok], s=28, color="#3d84b8", edgecolor="white", linewidth=0.4)
    ax3.plot(xx, intercept + slope * xx, color="#1f5f8b", linewidth=1.5)
    ax3.set_title("C. MYC dependency example", loc="left", weight="bold")
    ax3.set_xlabel("MTUS1 expression, log2(TPM + 1)")
    ax3.set_ylabel("MYC dependency (-Chronos score)")
    ax3.text(0.04, 0.95, f"n={ok.sum()}  r={r_value:.2f}, p={p_value:.1e}", transform=ax3.transAxes, va="top")

    ora = pd.DataFrame(
        {
            "pathway": [
                "MYC targets V1",
                "MYC targets V2",
                "Unfolded protein response",
                "G2M checkpoint",
                "E2F targets",
            ],
            "adjusted_p": [5.9e-4, 1.3e-3, 9.9e-3, 1.0e-2, 4.99e-2],
            "genes": [7, 4, 4, 5, 4],
        }
    )
    ora["score"] = -np.log10(ora["adjusted_p"])
    colors = ["#c94c4c", "#c94c4c", "#3f9f88", "#6b6fbd", "#6b6fbd"]
    ax4.barh(range(len(ora)), ora["score"], color=colors, alpha=0.9)
    ax4.set_yticks(range(len(ora)))
    ax4.set_yticklabels(ora["pathway"])
    ax4.invert_yaxis()
    ax4.axvline(-math.log10(0.05), color="#555555", linestyle="--", linewidth=0.8)
    ax4.set_xlabel("-log10(adjusted p-value)")
    ax4.set_title("D. Pathway convergence", loc="left", weight="bold")
    for i, row in ora.iterrows():
        ax4.text(row["score"] + 0.05, i, f"{row['genes']} genes", va="center", fontsize=7)

    for ax in [ax1, ax2, ax3, ax4]:
        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)
        ax.grid(axis="y", color="#eeeeee", linewidth=0.6)
        ax.set_axisbelow(True)

    save_figure(fig, "Figure2_DepMap_prioritisation")


def rounded_box(ax, xy, width, height, text, fc, ec="#333333", fontsize=9, weight="normal"):
    patch = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle="round,pad=0.02,rounding_size=0.025",
        linewidth=0.9,
        facecolor=fc,
        edgecolor=ec,
    )
    ax.add_patch(patch)
    ax.text(
        xy[0] + width / 2,
        xy[1] + height / 2,
        text.replace("\\n", "\n"),
        ha="center",
        va="center",
        fontsize=fontsize,
        weight=weight,
        linespacing=1.12,
    )
    return patch


def arrow(ax, start, end, color="#555555", style="-|>"):
    arr = FancyArrowPatch(start, end, arrowstyle=style, mutation_scale=9, linewidth=0.9, color=color)
    ax.add_patch(arr)


def plot_study_design() -> None:
    fig, ax = plt.subplots(figsize=(7.4, 3.2))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    boxes = [
        ((0.03, 0.60), 0.18, 0.20, "Public TNBC\\nRNA-seq cohorts", "#eef2f5"),
        ((0.28, 0.60), 0.18, 0.20, "MTUS1-low vs\\nMTUS1-high tertiles", "#eef2f5"),
        ((0.53, 0.60), 0.18, 0.20, "Differential\\nexpression + GSEA", "#eef2f5"),
        ((0.77, 0.60), 0.20, 0.20, "Reproducible\\ntumour programme", "#f8e9e4"),
        ((0.51, 0.24), 0.20, 0.17, "DepMap CRISPR\\nprioritisation", "#e8f1ee"),
        ((0.77, 0.24), 0.20, 0.17, "Published WEE1/\\nPKMYT1 support", "#f4f0e5"),
    ]
    for xy, w, h, text, color in boxes:
        rounded_box(ax, xy, w, h, text, color, fontsize=7.6, weight="bold" if "programme" in text else "normal")
    for x0, x1 in [(0.21, 0.28), (0.46, 0.53), (0.71, 0.77)]:
        arrow(ax, (x0, 0.70), (x1, 0.70))
    arrow(ax, (0.86, 0.60), (0.61, 0.41))
    arrow(ax, (0.86, 0.60), (0.87, 0.41))
    ax.text(0.02, 0.93, "Computational hypothesis-generation framework", fontsize=10.5, weight="bold")
    ax.text(0.56, 0.15, "Candidate modules for follow-up; not treatment validation", fontsize=7.5, color="#555555")
    save_figure(fig, "Figure0_Study_design")


def plot_working_model() -> None:
    fig, ax = plt.subplots(figsize=(7.4, 4.3))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    rounded_box(ax, (0.36, 0.43), 0.28, 0.17, "MTUS1-low\\nTNBC tumour profile", "#f8e9e4", weight="bold")
    modules = [
        ((0.05, 0.66), 0.25, 0.15, "MYC-linked\\ntranscriptional output", "#f3d7d2"),
        ((0.05, 0.22), 0.25, 0.15, "Proteostasis /\\nUPR", "#dceee7"),
        ((0.70, 0.66), 0.25, 0.15, "DNA repair /\\ncheckpoint stress", "#e2e4f2"),
        ((0.70, 0.22), 0.25, 0.15, "OXPHOS /\\nstress adaptation", "#f2e5d6"),
    ]
    for xy, w, h, text, color in modules:
        rounded_box(ax, xy, w, h, text, color, fontsize=9)
    arrow(ax, (0.36, 0.55), (0.30, 0.72))
    arrow(ax, (0.36, 0.45), (0.30, 0.30))
    arrow(ax, (0.64, 0.55), (0.70, 0.72))
    arrow(ax, (0.64, 0.45), (0.70, 0.30))
    rounded_box(ax, (0.31, 0.04), 0.38, 0.14, "External biological support:\\nWEE1/PKMYT1 findings in ATIP3-deficient models", "#f4f0e5", fontsize=8)
    arrow(ax, (0.50, 0.43), (0.50, 0.18), color="#777777")
    ax.text(0.03, 0.94, "Working model for experimental follow-up", fontsize=11, weight="bold")
    ax.text(0.03, 0.89, "The arrows denote hypothesis prioritisation, not established causal direction.", fontsize=8, color="#555555")
    save_figure(fig, "Figure3_Working_model")


def plot_proliferation_forest() -> None:
    stats = pd.DataFrame(
        {
            "cohort": ["SRP042620", "VUMC", "GSE192341", "GSE181466", "TCGA", "GSE202203", "SRP157974", "TNBC pooled"],
            "n": [28, 30, 34, 52, 78, 194, 298, 714],
            "pearson_r": [0.231, -0.141, -0.154, -0.404, -0.241, -0.536, -0.643, -0.611],
            "p": [0.236, 0.456, 0.386, 2.99e-3, 3.36e-2, 8.47e-16, 3.47e-36, 3.52e-74],
        }
    )
    stats.to_csv(OUT / "source_mtus1_proliferation_correlations.csv", index=False)
    fig, ax = plt.subplots(figsize=(4.4, 3.4))
    y = np.arange(len(stats))
    colors = np.where(stats["p"] < 0.05, "#3d84b8", "#b8c0ca")
    ax.scatter(stats["pearson_r"], y, s=stats["n"] / 2.5, color=colors, edgecolor="white", linewidth=0.5)
    ax.axvline(0, color="#777777", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{c} (n={n})" for c, n in zip(stats["cohort"], stats["n"])])
    ax.invert_yaxis()
    ax.set_xlabel("Pearson r: MTUS1 expression vs proliferation score")
    ax.set_title("MTUS1 expression and proliferation score", loc="left", weight="bold")
    ax.set_xlim(-0.75, 0.35)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.grid(axis="x", color="#eeeeee", linewidth=0.6)
    save_figure(fig, "FigureS1_Proliferation_correlation")


def main() -> None:
    nes, fdr = extract_tnbc_gsea()
    plot_study_design()
    plot_tnbc_core_dotplot(nes, fdr)
    expr, effects = read_depmap_data()
    corr = compute_depmap_correlations(expr, effects)
    plot_depmap_prioritisation(expr, effects, corr)
    plot_working_model()
    plot_proliferation_forest()
    print(f"Generated figures in {OUT}")


if __name__ == "__main__":
    main()
