from __future__ import annotations

import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import linregress


PROJECT = Path(__file__).resolve().parents[1]
DATA = PROJECT / "data" / "derived"
OUT = PROJECT / "figures"
OUT.mkdir(parents=True, exist_ok=True)


mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 8.3,
        "axes.titlesize": 9.5,
        "axes.labelsize": 8.2,
        "xtick.labelsize": 7.4,
        "ytick.labelsize": 7.4,
        "legend.fontsize": 7.4,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.7,
        "savefig.bbox": "tight",
    }
)


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.12,
        1.06,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=13,
        fontweight="bold",
        clip_on=False,
    )


def style_axes(ax: plt.Axes, axis: str = "y") -> None:
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.grid(axis=axis, color="#ececec", linewidth=0.6)
    ax.set_axisbelow(True)


def main() -> None:
    corr = pd.read_csv(DATA / "source_depmap_all_essentiality_correlations.csv")
    selected = corr[corr["selected_candidate"]].copy()
    myc = pd.read_csv(DATA / "depmap_myc_dependency_example.csv")

    fig = plt.figure(figsize=(7.5, 6.2))
    gs = fig.add_gridspec(
        2,
        2,
        height_ratios=[0.96, 1.0],
        width_ratios=[1.02, 1.00],
        hspace=0.43,
        wspace=0.54,
    )
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])

    bins = np.linspace(corr["Pearson Correlation"].min(), corr["Pearson Correlation"].max(), 42)
    ax1.hist(corr["Pearson Correlation"], bins=bins, color="#d9e0e8", edgecolor="white", linewidth=0.5)
    ax1.hist(selected["Pearson Correlation"], bins=bins, color="#c95f46", alpha=0.9)
    ax1.axvline(0, color="#777777", linewidth=0.8)
    ax1.axvline(0.4, color="#8f332c", linestyle="--", linewidth=0.9)
    shared_corr_xlim = (-0.8, 0.8)
    shared_corr_xticks = np.arange(-0.75, 0.76, 0.25)
    ax1.set_xlim(shared_corr_xlim)
    ax1.set_xticks(shared_corr_xticks)
    add_panel_label(ax1, "A.")
    ax1.set_xlabel("Pearson correlation r between\nChronos score and MTUS1 expression")
    ax1.set_ylabel("Genes")
    ax1.text(
        0.97,
        0.93,
        f"{len(selected)} candidates",
        transform=ax1.transAxes,
        ha="right",
        va="top",
        color="#8f332c",
    )
    style_axes(ax1)

    ax2.scatter(
        corr["Pearson Correlation"],
        corr["minus_log10_p"],
        s=11,
        color="#b8c2cf",
        edgecolor="none",
        alpha=0.72,
    )
    ax2.scatter(
        selected["Pearson Correlation"],
        selected["minus_log10_p"],
        s=20,
        color="#c95f46",
        edgecolor="white",
        linewidth=0.35,
        alpha=0.95,
        zorder=3,
    )
    ax2.axhline(-math.log10(0.05), color="#555555", linestyle="--", linewidth=0.8)
    ax2.axvline(0, color="#777777", linewidth=0.7)
    add_panel_label(ax2, "B.")
    ax2.set_xlabel("Pearson correlation r between\nChronos score and MTUS1 expression")
    ax2.set_ylabel("-log10(p)")
    ax2.set_xlim(shared_corr_xlim)
    ax2.set_xticks(shared_corr_xticks)
    ax2.set_ylim(0, max(3.8, corr["minus_log10_p"].max() + 0.45))

    label_positions = {
        "MYC": (0.73, 3.58, "left"),
        "DMAP1": (0.73, 3.28, "left"),
        "NOB1": (0.70, 2.96, "left"),
        "DDX49": (0.70, 2.78, "left"),
        "GGPS1": (0.70, 2.60, "left"),
        "PSMB6": (0.62, 2.08, "left"),
        "PSMB5": (0.62, 1.82, "left"),
        "HYOU1": (0.62, 1.58, "left"),
        "PRMT1": (0.62, 1.36, "left"),
        "POLR1C": (0.36, 1.70, "right"),
        "DPAGT1": (0.36, 1.44, "right"),
        "DKC1": (0.36, 1.18, "right"),
    }
    for gene, (x_text, y_text, ha) in label_positions.items():
        row = corr[corr["Gene"] == gene]
        if row.empty:
            continue
        r = row.iloc[0]
        ax2.annotate(
            gene,
            xy=(r["Pearson Correlation"], r["minus_log10_p"]),
            xytext=(x_text, y_text),
            textcoords="data",
            fontsize=6.5,
            ha=ha,
            va="center",
            bbox=None
            if gene == "DKC1"
            else {"boxstyle": "round,pad=0.16", "fc": "white", "ec": "none", "alpha": 0.94},
            arrowprops={"arrowstyle": "-", "color": "#777777", "lw": 0.55, "shrinkA": 0, "shrinkB": 2},
            clip_on=False,
            zorder=4,
        )
    style_axes(ax2)

    x = myc["MTUS1_expression"]
    y = myc["MYC_chronos"]
    slope, intercept, r_value, p_value, _ = linregress(x, y)
    xx = np.linspace(float(x.min()), float(x.max()), 100)
    ax3.scatter(x, y, s=28, color="#3d84b8", edgecolor="white", linewidth=0.45)
    ax3.plot(xx, intercept + slope * xx, color="#1f5f8b", linewidth=1.5)
    add_panel_label(ax3, "C.")
    ax3.set_xlabel("MTUS1 expression, log2(TPM + 1)")
    ax3.set_ylabel("MYC Chronos score")
    ax3.text(
        0.04,
        1.075,
        "MYC dependency",
        transform=ax3.transAxes,
        ha="left",
        va="top",
        fontsize=8.4,
        bbox={"boxstyle": "round,pad=0.18", "fc": "white", "ec": "none", "alpha": 0.88},
        clip_on=False,
    )
    ax3.text(
        0.04,
        1.015,
        f"n = {len(myc)}  r = {r_value:.2f}, p = {p_value:.1e}",
        transform=ax3.transAxes,
        ha="left",
        va="top",
        bbox={"boxstyle": "round,pad=0.18", "fc": "white", "ec": "none", "alpha": 0.88},
        clip_on=False,
    )
    style_axes(ax3)

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
            "group": ["MYC", "MYC", "Proteostasis", "Checkpoint", "Checkpoint"],
        }
    )
    colors = {"MYC": "#c94c4c", "Proteostasis": "#3f9f88", "Checkpoint": "#6b6fbd"}
    y_pos = np.arange(len(ora))
    ax4.hlines(y_pos, xmin=0.05, xmax=ora["adjusted_p"], color=[colors[g] for g in ora["group"]], linewidth=2.2, alpha=0.88)
    ax4.scatter(ora["adjusted_p"], y_pos, s=44, color=[colors[g] for g in ora["group"]], edgecolor="white", linewidth=0.5, zorder=3)
    ax4.set_yticks(range(len(ora)))
    ax4.set_yticklabels(["MYC targets V1", "MYC targets V2", "UPR", "G2M checkpoint", "E2F targets"])
    ax4.invert_yaxis()
    ax4.set_xscale("log")
    ax4.set_xlim(0.08, 4e-4)
    ax4.set_xticks([5e-2, 1e-2, 1e-3])
    ax4.set_xticklabels(["0.05", "0.01", "0.001"])
    ax4.axvline(0.05, color="#555555", linestyle="--", linewidth=0.8)
    ax4.set_xlabel("Adjusted p-value")
    add_panel_label(ax4, "D.")
    for i, row in ora.iterrows():
        ax4.text(row["adjusted_p"] * 0.82, i, f"{row['genes']} genes", va="center", fontsize=7.1)
    style_axes(ax4, axis="x")

    fig.savefig(OUT / "Figure2_DepMap_prioritisation.pdf")
    fig.savefig(OUT / "Figure2_DepMap_prioritisation.png", dpi=600)
    plt.close(fig)


if __name__ == "__main__":
    main()
