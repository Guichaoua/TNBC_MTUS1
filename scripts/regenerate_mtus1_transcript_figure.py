from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
DATA = PROJECT / "data" / "derived" / "mtus1_transcript_summary.csv"
OUT = PROJECT / "figures"
OUT.mkdir(parents=True, exist_ok=True)


mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 8.4,
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


ORDER = ["ATIP1", "ATIP3a_total", "ATIP3b_total", "ATIP4"]
COHORTS = ["TCGA BRCA TNBC", "TNBC cell lines"]
COLORS = {
    "ATIP1": "#4c78a8",
    "ATIP3a_total": "#d76545",
    "ATIP3b_total": "#54a24b",
    "ATIP4": "#8e6bbd",
}


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.12, 1.06, label, transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")


def style_axes(ax: plt.Axes, axis: str = "y") -> None:
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.grid(axis=axis, color="#ececec", linewidth=0.6)
    ax.set_axisbelow(True)


def main() -> None:
    df = pd.read_csv(DATA)
    df["transcript_group"] = pd.Categorical(df["transcript_group"], ORDER, ordered=True)
    df["cohort"] = pd.Categorical(df["cohort"], COHORTS, ordered=True)
    df = df.sort_values(["cohort", "transcript_group"])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.45, 3.25), gridspec_kw={"width_ratios": [1.15, 1.0]})

    width = 0.18
    x = np.arange(len(COHORTS))
    for i, transcript in enumerate(ORDER):
        sub = df[df["transcript_group"] == transcript].set_index("cohort").loc[COHORTS]
        positions = x + (i - 1.5) * width
        ax1.bar(positions, sub["mean_tpm"], width=width, color=COLORS[transcript], label=transcript)
        for xpos, value in zip(positions, sub["mean_tpm"]):
            if value >= 0.15:
                ax1.text(xpos, value + 0.10, f"{value:.1f}", ha="center", va="bottom", fontsize=6.7)

    add_panel_label(ax1, "A")
    ax1.set_title("Mean transcript expression", loc="left", fontweight="bold")
    ax1.set_ylabel("Mean TPM")
    ax1.set_xticks(x)
    ax1.set_xticklabels(["TCGA BRCA\nTNBC tumours", "TNBC\ncell lines"])
    ax1.set_ylim(0, max(df["mean_tpm"]) * 1.26)
    style_axes(ax1)

    left = np.zeros(len(COHORTS))
    for transcript in ORDER:
        sub = df[df["transcript_group"] == transcript].set_index("cohort").loc[COHORTS]
        values = sub["percent_total"].to_numpy()
        ax2.barh(COHORTS, values, left=left, color=COLORS[transcript], label=transcript)
        for j, value in enumerate(values):
            if value >= 8:
                ax2.text(left[j] + value / 2, j, f"{value:.1f}%", ha="center", va="center", fontsize=7, color="white")
            elif value >= 2:
                ax2.text(left[j] + value + 1.2, j, f"{value:.1f}%", ha="left", va="center", fontsize=6.6, color="#333333")
        left += values

    add_panel_label(ax2, "B")
    ax2.set_title("Retained coding MTUS1 signal", loc="left", fontweight="bold")
    ax2.set_xlabel("Contribution to retained signal (%)")
    ax2.set_xlim(0, 100)
    ax2.set_yticks(range(len(COHORTS)))
    ax2.set_yticklabels(["TCGA BRCA\nTNBC tumours", "TNBC\ncell lines"])
    style_axes(ax2, axis="x")

    handles, labels = ax1.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.53, -0.06))
    fig.subplots_adjust(bottom=0.24, left=0.08, right=0.98, wspace=0.34)

    fig.savefig(OUT / "MTUS1_transcripts_article.pdf")
    fig.savefig(OUT / "MTUS1_transcripts_article.png", dpi=600)
    plt.close(fig)


if __name__ == "__main__":
    main()
