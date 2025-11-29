import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def make_plots(out_dir: str):
    """
    Generate clinician-friendly PNG plots from the TSV outputs.
    Includes:
      - SV type distribution
      - Inheritance categories
      - Clinically relevant annotation overlaps
      - Pathogenicity tier distribution
      - De novo SV overview
      - Per-chromosome SV counts by type
      - SV length distribution (log10)
    """

    def safe_read_tsv(path, **kwargs):
        if not os.path.exists(path):
            print(f"[WARN] File not found for plotting: {path}")
            return None
        return pd.read_csv(path, sep="\t", **kwargs)

    # ----------------------------------------------------------
    # FIGURE 1: SVTYPE composition (DEL / INS / DUP)
    # ----------------------------------------------------------
    svtype_path = os.path.join(out_dir, "svtype_counts.tsv")
    sv_counts = safe_read_tsv(svtype_path)
    if sv_counts is not None and {"SVTYPE", "COUNT"}.issubset(sv_counts.columns):
        sv_counts = sv_counts.sort_values("COUNT", ascending=False)

        plt.figure(figsize=(5, 4))
        plt.bar(sv_counts["SVTYPE"], sv_counts["COUNT"])
        for x, y in zip(sv_counts["SVTYPE"], sv_counts["COUNT"]):
            plt.text(x, y, str(y), ha="center", va="bottom", fontsize=8)
        plt.title("Structural Variant Types")
        plt.xlabel("SV type")
        plt.ylabel("Number of SVs")
        plt.tight_layout()
        fig1_path = os.path.join(out_dir, "fig1_svtype_counts.png")
        plt.savefig(fig1_path, dpi=300)
        plt.close()
        print(f"[INFO] Saved {fig1_path}")

    # ----------------------------------------------------------
    # FIGURE 2: Inheritance categories (child-centric)
    # ----------------------------------------------------------
    inh_path = os.path.join(out_dir, "inheritance_counts.tsv")
    inh_counts = safe_read_tsv(inh_path)
    if inh_counts is not None and {"INHERITANCE", "COUNT"}.issubset(inh_counts.columns):
        order = [
            "de_novo",
            "inherited_p1",
            "inherited_p2",
            "inherited_both",
            "reference_or_absent_child",
            "unknown",
        ]
        cats_in_data = [c for c in order if c in inh_counts["INHERITANCE"].tolist()]
        others = [c for c in inh_counts["INHERITANCE"].tolist() if c not in cats_in_data]
        ordered = cats_in_data + others
        inh_counts["INHERITANCE"] = pd.Categorical(
            inh_counts["INHERITANCE"], categories=ordered, ordered=True
        )
        inh_counts = inh_counts.sort_values("INHERITANCE")

        plt.figure(figsize=(7, 4))
        plt.bar(inh_counts["INHERITANCE"].astype(str), inh_counts["COUNT"])
        for x, y in zip(inh_counts["INHERITANCE"].astype(str), inh_counts["COUNT"]):
            plt.text(x, y, str(y), ha="center", va="bottom", fontsize=8)
        plt.title("Inheritance Pattern of SVs in the Child")
        plt.xlabel("Inheritance category")
        plt.ylabel("Number of SVs")
        plt.xticks(rotation=30, ha="right")
        plt.tight_layout()
        fig2_path = os.path.join(out_dir, "fig2_inheritance_counts.png")
        plt.savefig(fig2_path, dpi=300)
        plt.close()
        print(f"[INFO] Saved {fig2_path}")

    # ----------------------------------------------------------
    # FIGURE 3: Annotation overlaps (Exons, HI, TS, Recurrent CNV, ClinVar)
    # uses *_counts.tsv files (bedtools -c outputs)
    # ----------------------------------------------------------
    annot_specs = [
        ("Exon-overlapping SVs", "hg38_exons_counts.tsv"),
        ("ClinGen HI genes", "clingen_hi_counts.tsv"),
        ("ClinGen TS genes", "clingen_ts_counts.tsv"),
        ("Recurrent CNV regions", "clingen_recurrent_cnv_counts.tsv"),
        ("ClinVar pathogenic SVs", "clinvar_sv_counts.tsv"),
    ]

    annot_rows = []
    for label, fname in annot_specs:
        path = os.path.join(out_dir, fname)
        if not os.path.exists(path):
            print(f"[WARN] Missing annotation counts file for {label}: {path}")
            continue
        df_raw = pd.read_csv(path, sep="\t", header=None)
        if df_raw.shape[1] < 5:
            print(f"[WARN] Unexpected format in {path}; skipping.")
            continue
        n_overlap = (df_raw.iloc[:, 4] > 0).sum()
        annot_rows.append({"category": label, "n_overlap": int(n_overlap)})

    if annot_rows:
        annot_df = pd.DataFrame(annot_rows).sort_values("n_overlap", ascending=True)

        plt.figure(figsize=(7, 4.5))
        plt.barh(annot_df["category"], annot_df["n_overlap"])
        for y, val in zip(annot_df["category"], annot_df["n_overlap"]):
            plt.text(val, y, str(val), va="center", ha="left", fontsize=8)
        plt.title("Clinically Relevant Annotations")
        plt.xlabel("Number of SVs with ≥1 overlap")
        plt.ylabel("")
        plt.tight_layout()
        fig3_path = os.path.join(out_dir, "fig3_annotation_overlap.png")
        plt.savefig(fig3_path, dpi=300)
        plt.close()
        print(f"[INFO] Saved {fig3_path}")

    # ----------------------------------------------------------
    # FIGURE 4: Pathogenicity tiers of clinically prioritised SVs
    # ----------------------------------------------------------
    voi_path = os.path.join(out_dir, "clinically_prioritised_SVs.tsv")
    voi = safe_read_tsv(voi_path)
    if voi is not None and "PATHOGENICITY_TIER" in voi.columns:
        tier_counts = (
            voi["PATHOGENICITY_TIER"]
            .value_counts()
            .rename_axis("PATHOGENICITY_TIER")
            .reset_index(name="COUNT")
            .sort_values("COUNT", ascending=False)
        )

        plt.figure(figsize=(7, 4))
        plt.bar(tier_counts["PATHOGENICITY_TIER"], tier_counts["COUNT"])
        for x, y in zip(tier_counts["PATHOGENICITY_TIER"], tier_counts["COUNT"]):
            plt.text(x, y, str(y), ha="center", va="bottom", fontsize=8)
        plt.title("Clinically Prioritised SVs by Pathogenicity Tier")
        plt.xlabel("Pathogenicity tier")
        plt.ylabel("Number of SVs")
        plt.xticks(rotation=30, ha="right")
        plt.tight_layout()
        fig4_path = os.path.join(out_dir, "fig4_pathogenicity_tiers.png")
        plt.savefig(fig4_path, dpi=300)
        plt.close()
        print(f"[INFO] Saved {fig4_path}")

    # ----------------------------------------------------------
    # FIGURE 5: Per-chromosome SV counts by SVTYPE
    # ----------------------------------------------------------
    trio_table_path = os.path.join(out_dir, "trio_table.tsv")
    trio = safe_read_tsv(trio_table_path)
    if trio is not None and {"CHROM", "SVTYPE"}.issubset(trio.columns):
        chrom_sv = (
            trio.groupby(["CHROM", "SVTYPE"])
            .size()
            .rename("COUNT")
            .reset_index()
        )

        # consistent chromosome order
        def chrom_order(c):
            c = str(c)
            if c.startswith("chr"):
                c2 = c[3:]
            else:
                c2 = c
            try:
                return (0, int(c2))
            except ValueError:
                # X, Y, etc.
                if c2 == "X":
                    return (1, 23)
                if c2 == "Y":
                    return (1, 24)
                return (2, 999)
        chroms = sorted(chrom_sv["CHROM"].unique(), key=chrom_order)

        svtypes = sorted(chrom_sv["SVTYPE"].unique())
        x = np.arange(len(chroms))
        width = 0.8 / max(1, len(svtypes))

        plt.figure(figsize=(8, 4))
        for i, svt in enumerate(svtypes):
            sub = chrom_sv[chrom_sv["SVTYPE"] == svt]
            counts = [sub[sub["CHROM"] == c]["COUNT"].sum() for c in chroms]
            plt.bar(x + i * width, counts, width=width, label=str(svt))

        plt.xticks(x + width * (len(svtypes) - 1) / 2, chroms, rotation=45, ha="right")
        plt.ylabel("Number of SVs")
        plt.xlabel("Chromosome")
        plt.title("SV Counts per Chromosome by Type")
        plt.legend(title="SVTYPE", fontsize=8)
        plt.tight_layout()
        fig6_path = os.path.join(out_dir, "fig5_per_chrom_svtype.png")
        plt.savefig(fig6_path, dpi=300)
        plt.close()
        print(f"[INFO] Saved {fig6_path}")

    # ----------------------------------------------------------
    # FIGURE 7: SV length distribution (log10), for child SVs
    # ----------------------------------------------------------
    if trio is not None and "SVLEN_BP" in trio.columns:
        plt.figure(figsize=(6, 4))
        lengths = trio["SVLEN_BP"].replace(0, np.nan).dropna().abs()
        plt.hist(np.log10(lengths), bins=40)
        plt.title("SV Length Distribution (log10 scale)")
        plt.xlabel("log10(SV length in bp)")
        plt.ylabel("Count")
        plt.tight_layout()
        fig7_path = os.path.join(out_dir, "fig7_sv_length_log10.png")
        plt.savefig(fig7_path, dpi=300)
        plt.close()
        print(f"[INFO] Saved {fig7_path}")
