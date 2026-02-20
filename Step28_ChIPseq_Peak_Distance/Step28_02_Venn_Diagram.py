#!/usr/bin/env python3
"""
Python 3.12+ script to:
1) Read two BED files (chr, start, end).
2) Compute an overlap-based Venn partition (A-only, B-only, shared).
3) Render an area-proportional Venn diagram WITHOUT count numbers.
4) Save it to:
   ../../public/Multiomics/Step28_ChIPseq_Peak_Distance/Peak_Distance_Venn.png

Dependencies:
  pip install pyranges pandas matplotlib matplotlib-venn
"""

from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import pyranges as pr


def read_bed_as_pyranges(bed_path: Path) -> pr.PyRanges:
    if not bed_path.exists():
        raise FileNotFoundError(f"Missing BED file: {bed_path}")
    # pr.read_bed reads common BED formats into a PyRanges object.
    return pr.read_bed(str(bed_path))


def venn_counts_by_cluster(a: pr.PyRanges, b: pr.PyRanges, label_a: str, label_b: str) -> tuple[int, int, int]:
    """
    Produces symmetric A-only / B-only / shared counts by:
    - tagging intervals by source
    - concatenating
    - clustering by genomic overlap
    - counting clusters that contain only A, only B, or both

    Note: Using slack=-1 makes "bookended" intervals (End == Start) NOT count as overlap,
    which typically matches BED's half-open semantics.
    """
    df_a = a.df.copy()
    df_b = b.df.copy()
    df_a["Source"] = label_a
    df_b["Source"] = label_b

    combined = pr.PyRanges(pd.concat([df_a, df_b], ignore_index=True))

    # Cluster overlapping intervals across both sets.
    clustered = combined.cluster(slack=-1)
    cdf = clustered.df

    # For each cluster, determine which sources appear in it.
    sources_by_cluster = (
        cdf.groupby("Cluster")["Source"]
        .apply(lambda s: frozenset(pd.unique(s)))
    )

    a_only = int((sources_by_cluster == frozenset([label_a])).sum())
    b_only = int((sources_by_cluster == frozenset([label_b])).sum())
    shared = int((sources_by_cluster == frozenset([label_a, label_b])).sum())

    return a_only, b_only, shared


def main() -> int:
    bed_a = Path("../../public/ChIPseq/SHK_2024_06_11_MCF7_DBN021_1_0_ug_mL_Rep2 Filtered_trim.bed")
    bed_b = Path("../../public/ChIPseq/SHK_2024_06_11_MCF7_DBN025_1_0_ug_mL_Rep2 Filtered_trim.bed")
    out_png = Path("../../public/Multiomics/Step28_ChIPseq_Peak_Distance/Peak_Distance_Venn.png")
    out_png.parent.mkdir(parents=True, exist_ok=True)

    label_a = "1"
    label_b = "2"

    # Read BED files
    try:
        a = read_bed_as_pyranges(bed_a)
        b = read_bed_as_pyranges(bed_b)
    except FileNotFoundError as e:
        print(e, file=sys.stderr)
        return 1

    # Calculate overlaps
    a_only, b_only, shared = venn_counts_by_cluster(a, b, label_a=label_a, label_b=label_b)

    # Setup plot
    fig, ax = plt.subplots(figsize=(6.5, 6.0))
    
    # Draw Venn diagram
    v = venn2(subsets=(a_only, b_only, shared), set_labels=(label_a, label_b), ax=ax)
    
    # Iterate through the subset labels (the count numbers) and hide them
    if v.subset_labels:
        for label in v.subset_labels:
            if label:
                label.set_visible(False)

    ax.set_title("ChIP-seq peaks overlap (clustered by genomic overlap)")

    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved Venn diagram (without counts) to: {out_png}")
    print(f"Calculated counts - A-only: {a_only}, B-only: {b_only}, Shared: {shared}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())