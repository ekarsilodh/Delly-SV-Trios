#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Author: Ekarsi Lodh

"""
VCF Trio Analyzer (AWK + Python)
--------------------------------

Pipeline:
1. AWK pass:
    - Reads VCF
    - Extracts CHROM..SVLEN, and GT for the LAST THREE samples
    - Fast parsing

2. Python:
    - Infers trio roles (child / parent1 / parent2)
    - Determines ORIGINAL COLUMN numbers (exact like $16, $17, $18 etc)
    - Runs AWK de-novo finder 
    - Summarizes:

       a) De-novo SV calls
       b) All trio calls (counts SVTYPE, etc.)

    - bedtools annotation step
       * gene overlap
       * exon overlap
       * ClinGen HI / TS / recurrent CNVs

Outputs:
    - trio_table.tsv              : Parsed AWK table with GTs and SV annotations.
    - svtype_counts.tsv           : Global count of SVTYPEs.
    - de_novo_calls.txt / .tsv    : De novo candidate SVs (0/0 parents, child non-ref).
    - clinically_prioritised_SVs.tsv : Child's SVs passing size/region filters.
    - sv_for_bedtools.bed         : BED intervals used for all bedtools annotations.
    - bedtools_hg38_genes_raw.tsv : Raw SV × gene overlaps.
    - hg38_exons_raw.tsv / _counts.tsv : Raw + counts for exon overlaps.
    - clingen_hi_raw.tsv / _counts.tsv : Raw + counts for ClinGen HI genes.
    - clingen_ts_raw.tsv / _counts.tsv : Raw + counts for ClinGen TS genes.
    - clingen_recurrent_cnv_raw.tsv / _counts.tsv : Raw + counts for recurrent CNVs.
    - report.txt                  : Compact summary.
    - summary_full.txt            : Full overview of methods and findings.
"""

import argparse
import subprocess
import tempfile
import pandas as pd
import os
import sys
from itertools import permutations
from typing import Optional, Dict

#custom module for plotting the results
from plot import make_plots


############################################
# CONSTANTS / PATHS FOR BEDTOOLS ANNOTATION
############################################

HG38_GENES_BED = "../databases/hg38_genes.bed"
HG38_EXONS_BED = "../databases/hg38_exons.bed"

CLINGEN_HI_GENES_BED = "../databases/ClinGen_haploinsufficiency_gene_GRCh38.bed"
CLINGEN_TS_GENES_BED = "../databases/ClinGen_triplosensitivity_gene_GRCh38.bed"
CLINGEN_RECURRENT_CNV_BED = "../databases/ClinGen recurrent CNV-hg38.bed"

CLINVAR_SV_BED = "clinvar_SV_clean.bed"  # Pathogenic / Likely pathogenic SVs (GRCh38) only

# 1. ARGUMENT PARSING

def parse_args():
    """
    Parse command-line arguments.

    Returns a Namespace with:
      vcf: path to input VCF
      out: output directory (default: "out")
    """
    p = argparse.ArgumentParser()
    p.add_argument("--vcf", required=True)
    p.add_argument("--out", default="../out")
    return p.parse_args()

def get_vcf_samples(vcf_path: str):
    """
    Return the list of sample names from the VCF header (columns 10+ on the #CHROM line).
    """
    header_samples = []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#CHROM"):
                parts = line.strip().split("\t")
                header_samples = parts[9:]
                break
    if not header_samples:
        raise RuntimeError("Could not find any samples in VCF header (#CHROM line).")
    return header_samples


def choose_trio_samples(vcf_path: str):
    """
    Interactively show the sample names in the VCF and ask the user
    to type the three sample names they want to analyse.
    Returns a list [sample1, sample2, sample3].
    """
    all_samples = get_vcf_samples(vcf_path)

    print("\n===== SAMPLES IN VCF =====")
    for i, s in enumerate(all_samples, start=1):
        print(f"{i}. {s}")

    print("\nEnter the names of the 3 samples you want to analyse,")
    print("exactly as they appear above.")

    chosen = []
    while len(chosen) < 3:
        name = input(f"Sample {len(chosen)+1} name: ").strip()
        if name not in all_samples:
            print(f"  → '{name}' is not a valid sample name in this VCF. Please try again.")
            continue
        if name in chosen:
            print("  → You already selected that sample; please choose a different one.")
            continue
        chosen.append(name)

    return chosen


# 2. AWK SCRIPT FOR FIRST PASS
# awk: pick last 3 samples as trio, pull out GT + SV info

AWK_PARSE = r'''
function extract_gt(sampleField,   a,b) {
    n = split(sampleField, a, ":");
    # return genotype only
    split(a[1], b, ":");
    return a[1];
}

BEGIN { OFS="\t"; }

/^##/ { next }

/^#CHROM/ {
    for (i=1; i<=NF; i++) H[i]=$i;
    n=NF;

    # Find the three sample columns matching the requested names.
    # The variables sample1, sample2 and sample3 are passed in from Python via -v.
    s1 = s2 = s3 = -1;
    for (i=10; i<=NF; i++) {
        if (H[i] == sample1)      s1 = i;
        else if (H[i] == sample2) s2 = i;
        else if (H[i] == sample3) s3 = i;
    }

    # Fallback: if any of the requested names were not found,
    # default to using the last three sample columns.
    if (s1 < 0 || s2 < 0 || s3 < 0) {
        s1 = n-2; s2 = n-1; s3 = n;
    }

    print "CHROM","POS","REF","ALT","QUAL","FILTER","INFO",
          H[s1]"_GT",H[s2]"_GT",H[s3]"_GT", "SVTYPE", "END", "SVLEN";
    next;
}


/^[^#]/ && $7 == "PASS" {
    # parse INFO
    SVTYPE=""; ENDV=""; SVLEN="";
    split($8, infoA, ";");
    for (i in infoA){
        split(infoA[i], kv, "=");
        if (kv[1]=="SVTYPE") SVTYPE=kv[2];
        else if (kv[1]=="END") ENDV=kv[2];
        else if (kv[1]=="SVLEN") SVLEN=kv[2];
    }

    gt1 = extract_gt($(s1));
    gt2 = extract_gt($(s2));
    gt3 = extract_gt($(s3));

    print $1,$2,$4,$5,$6,$7,$8,gt1,gt2,gt3,SVTYPE,ENDV,SVLEN;
}
'''


# 3. GENOTYPE CODING

def gt_to_code(gt: Optional[str]) -> Optional[int]:
    """Convert GT to simple 0/1/2 code."""
    if gt in ("0/0", "0|0", "0"):
        return 0
    if gt in ("0/1", "1/0", "0|1", "1|0"):
        return 1
    if gt in ("1/1", "1|1", "1"):
        return 2
    return None


# 4. TRIO ROLE INFERENCE by MENDELIAN VIOLATION

def mendelian_violation(child, p1, p2):
    """
    - Parent1 genotype = 0/0  → encoded as 0
    - Parent2 genotype = 0/0  → encoded as 0
    - Child genotype  != 0/0  → encoded as anything non-zero

    Returns True iff this pattern is observed, False otherwise.
    """
    # ignore if any genotype is missing
    if child is None or p1 is None or p2 is None:
        return False

    # both parents 0/0, child not 0/0
    return (p1 == 0 and p2 == 0 and child != 0)


def infer_trio(df, samples):
    """
    Given three samples, try all child/parent assignments and
    pick the one with the fewest Mendelian violations.
    """
    best = None
    roles = None
    trio_stats = []

    for child, p1, p2 in permutations(samples, 3):
        vio = 0
        for _, row in df.iterrows():
            c = gt_to_code(row[f"{child}_GT"])
            m = gt_to_code(row[f"{p1}_GT"])
            f = gt_to_code(row[f"{p2}_GT"])
            if mendelian_violation(c, m, f):
                vio += 1

        trio_stats.append({
            "child": child,
            "parent1": p1,
            "parent2": p2,
            "violations": vio,
        })

        if best is None or vio < best:
            best = vio
            roles = {"child": child, "parent1": p1, "parent2": p2}
    
    # sanity check — if even the best assignment has > 4 violations,
    # this trio probably isn't a valid family.
    if best is not None and best > 4:
        print("\n[ERROR] Minimum Mendelian-violation count across all trio assignments "
              f"is {best} (> 4).")
        print("[ERROR] This trio does not look like a valid family. Aborting analysis.")
        sys.exit(1)
        
    return roles, trio_stats


# 5. SEX INFERENCE FROM chrX

def infer_sex_from_chrX(
    df: pd.DataFrame,
    samples,
    het_threshold: float = 0.1,
    min_sites: int = 5,
):
    """
    Sex inference based on heterozygosity on chrX:

      - For each sample:
          H = (# het genotypes on X) / (# non-missing genotypes on X)
      - If H < het_threshold → "M"
      - Else → "F"
    Returns:
      sex_calls: {sample: "M"/"F"/"unknown"}
      het_rates: {sample: float or None}
    """

    sex_calls: Dict[str, str] = {s: "unknown" for s in samples}
    het_rates: Dict[str, Optional[float]] = {s: None for s in samples}


    chrom_series = df["CHROM"].astype(str)
    chrom_nochr = chrom_series.str.replace("^chr", "", case=False, regex=True)
    mask_x = chrom_nochr == "X"
    df_x = df[mask_x]

    if df_x.empty:
        return sex_calls, het_rates

    for s in samples:
        col = f"{s}_GT"
        if col not in df_x.columns:
            continue

        gts = df_x[col].astype(str)
        non_missing = ~gts.isin(["", ".", "./.", ".|."])
        n_total = int(non_missing.sum())
        if n_total < min_sites:
            continue

        gts_nm = gts[non_missing]
        het = int(gts_nm.isin(["0/1", "1/0", "0|1", "1|0"]).sum())
        H = het / n_total if n_total > 0 else 0.0

        het_rates[s] = H
        if H < het_threshold:
            sex_calls[s] = "M"
        else:
            sex_calls[s] = "F"

    return sex_calls, het_rates

# 6. DE NOVO DETECTION VIA AWK

def run_awk_denovo(vcf, outpath, outpath1, parent1_col, parent2_col, child_col):
    """
    Use awk on the original VCF to pull out obvious de novo calls.

    Parents must be 0/0, child must be non-ref and not missing.
    Writes two identical tab-separated files (outpath, outpath1).
    """

    awk_expr = (
        'BEGIN{OFS="\\t"; '
        'print "CHROM","START_POS","END_POS","ID","REF","ALT","SV_TYPE","LENGTH_bp","Paired_end_PE","Split_end_SR",'
        f'"Mother_GT","Father_GT","Child_GT"; }} '
        f'$0!~/^#/ && ${parent1_col} ~ /^0\\/0/ && '
        f'${parent2_col} ~ /^0\\/0/ && '
        f'${child_col} !~ /^0\\/0/ && '
        f'${child_col} !~ /^\\.\\/\\./ '
        '{ '
            'pos = $2; '

            # INFO PARSE
            'n = split($8, element, ";"); '
            'for (i = 1; i <= n; i++) { '
                'split(element[i], kv, "="); '
                'key = kv[1]; val = kv[2]; '
                'if (key == "SVTYPE") svtype = val; '
                'else if (key == "END") endpos = val + 0; '
                'else if (key == "PE") pe = val + 0; '
                'else if (key == "SR") sr = val + 0; '
            '} '

            # LENGTH
            'length_bp = endpos - pos; '

            # GENOTYPES
            f'split(${parent1_col}, a, ":"); g1 = a[1]; '
            f'split(${parent2_col}, b, ":"); g2 = b[1]; '
            f'split(${child_col}, c, ":"); g3 = c[1]; '

            # OUTPUT
            'print $1,pos,endpos,$3,$4,$5,svtype,length_bp,pe,sr,g1,g2,g3; '
        '}'
    )

    cmd = ["awk", awk_expr, vcf]

    # Run awk and capture to file 1
    with open(outpath, "w") as outfh:
        subprocess.run(cmd, check=True, stdout=outfh)

    # Run awk and capture to file 2 (same content, different path)
    with open(outpath1, "w") as outfh:
        subprocess.run(cmd, check=True, stdout=outfh)

        
# 7. SV LENGTH, SIZE PRIORITY

def add_sv_length_and_size_category(df: pd.DataFrame) -> pd.DataFrame:
    """
    Adds:
      - POS_INT, END_INT, SVLEN_INT (numeric)
      - SVLEN_BP      : final length used (abs SVLEN or END - POS)
      - SIZE_PRIORITY : 'high', 'moderate', 'low', or 'unknown'
    This is applied to ALL SVs, but is most meaningful for DEL/DUP.
    """
    # numeric POS
    df["POS_INT"] = pd.to_numeric(df["POS"], errors="coerce")

    # numeric END (from parsed TSV)
    if "END" in df.columns:
        df["END_INT"] = pd.to_numeric(df["END"], errors="coerce")
    else:
        df["END_INT"] = pd.NA

    # numeric SVLEN
    if "SVLEN" in df.columns:
        df["SVLEN_INT"] = pd.to_numeric(df["SVLEN"], errors="coerce")
    else:
        df["SVLEN_INT"] = pd.NA

    # prefer SVLEN, otherwise END - POS
    def compute_len(row):
        if pd.notna(row.get("SVLEN_INT")):
            try:
                return abs(int(row["SVLEN_INT"]))
            except Exception:
                return pd.NA
        if pd.notna(row.get("END_INT")) and pd.notna(row.get("POS_INT")):
            try:
                return max(0, int(row["END_INT"]) - int(row["POS_INT"]))
            except Exception:
                return pd.NA
        return pd.NA

    df["SVLEN_BP"] = df.apply(compute_len, axis=1)

    def size_priority(length_bp):
        """
        Simple thresholds:
          >100 kb       → high
          20–100 kb     → moderate
          <20 kb        → low
        """
        if pd.isna(length_bp):
            return "unknown"
        length_bp = int(length_bp)
        if length_bp > 100_000:
            return "high"
        elif 20_000 <= length_bp <= 100_000:
            return "moderate"
        else:
            return "low"

    df["SIZE_PRIORITY"] = df["SVLEN_BP"].apply(size_priority)
    return df


# 8. INHERITANCE LABEL

def annotate_inheritance(df: pd.DataFrame, child: str, p1: str, p2: str) -> pd.DataFrame:
    """
    Label inheritance for ALL variants based on trio genotypes:

      - de_novo
      - inherited_p1
      - inherited_p2
      - inherited_both
      - reference_or_absent_child
      - unknown

    This does NOT filter anything; it just adds an INHERITANCE column.
    """
    child_col = f"{child}_GT"
    p1_col = f"{p1}_GT"
    p2_col = f"{p2}_GT"

    def norm(gt):
        if not isinstance(gt, str):
            return ""
        return gt.strip()

    def is_ref(gt):
        return gt in ("0/0", "0|0", "0")

    def is_alt(gt):
        # any non-ref diploid genotype
        return gt in ("0/1", "1/0", "0|1", "1|0", "1/1", "1|1", "1")

    def is_missing(gt):
        return gt in ("", ".", "./.", ".|.")

    def classify(row):
        c = norm(row.get(child_col, ""))
        m = norm(row.get(p1_col, ""))
        f = norm(row.get(p2_col, ""))

        if is_missing(c) or (is_missing(m) and is_missing(f)):
            return "unknown"

        if is_ref(c) or not is_alt(c):
            return "reference_or_absent_child"

        m_alt = is_alt(m)
        f_alt = is_alt(f)

        if is_ref(m) and is_ref(f):
            return "de_novo"
        if m_alt and not f_alt:
            return "inherited_p1"
        if f_alt and not m_alt:
            return "inherited_p2"
        if m_alt and f_alt:
            return "inherited_both"
        return "unknown"

    df["INHERITANCE"] = df.apply(classify, axis=1)
    return df


# 9. Annotation Using bedtools

def bedtools_annotate_all_variants(df: pd.DataFrame, out_dir: str) -> pd.DataFrame:
    """
    Use bedtools intersect to annotate every SV with overlaps to:
      - hg38 gene bodies (HG38_GENES_BED) -> GENE_OVERLAP_COUNT and GENE_LIST
      - hg38 exons (HG38_EXONS_BED) -> EXON_OVERLAP_COUNT
      - ClinGen HI genes (CLINGEN_HI_GENES_BED) -> CLINGEN_HI_OVERLAP_COUNT
      - ClinGen TS genes (CLINGEN_TS_GENES_BED) -> CLINGEN_TS_OVERLAP_COUNT
      - ClinGen recurrent CNVs (CLINGEN_RECURRENT_CNV_BED) -> CLINGEN_REC_CNV_OVERLAP_COUNT

    All coordinates are assumed to be hg38/GRCh38, matching the input VCF.
    Raw bedtools outputs are written under out_dir for inspection.
    """

    df = df.copy()

    # Ensure numeric positions/ends
    if "POS_INT" not in df.columns:
        df["POS_INT"] = pd.to_numeric(df["POS"], errors="coerce")

    if "END_INT" not in df.columns:
        if "END" in df.columns:
            df["END_INT"] = pd.to_numeric(df["END"], errors="coerce")
        else:
            df["END_INT"] = df["POS_INT"]

    # Create a unique row id to map back bedtools results
    df["ROW_ID"] = range(len(df))

    bed_input = os.path.join(out_dir, "sv_for_bedtools.bed")
    with open(bed_input, "w") as bed_out:
        for _, row in df.iterrows():
            if pd.isna(row["POS_INT"]):
                continue
            chrom = str(row["CHROM"])
            try:
                pos = int(row["POS_INT"])
            except Exception:
                continue
            end = row["END_INT"]
            if pd.isna(end):
                end = pos
            else:
                try:
                    end = int(end)
                except Exception:
                    end = pos
            # BED: 0-based start, 1-based end
            start0 = max(pos - 1, 0)
            bed_out.write(f"{chrom}\t{start0}\t{end}\t{int(row['ROW_ID'])}\n")

    def annotate_bed(bed_path: str, count_col: str, prefix: str):
        """
        Run bedtools -c and store counts in df[count_col].
        """
        if not bed_path:
            df[count_col] = 0
            return

        if not os.path.exists(bed_path):
            print(f"[WARN] BED file not found: {bed_path}; skipping {count_col}.")
            df[count_col] = 0
            return

        try:
            # Raw overlaps (-wa -wb)
            raw_file = os.path.join(out_dir, f"{prefix}_raw.tsv")
            with open(raw_file, "w") as raw_fh:
                subprocess.run(
                    ["bedtools", "intersect", "-a", bed_input, "-b", bed_path, "-wa", "-wb"],
                    stdout=raw_fh,
                    check=True,
                )

            # Counts (-c)
            counts_file = os.path.join(out_dir, f"{prefix}_counts.tsv")
            with open(counts_file, "w") as cf:
                subprocess.run(
                    ["bedtools", "intersect", "-a", bed_input, "-b", bed_path, "-c"],
                    stdout=cf,
                    check=True,
                )

            counts = {}
            with open(counts_file, "r") as cf:
                for line in cf:
                    parts = line.rstrip().split("\t")
                    if len(parts) < 5:
                        continue
                    try:
                        row_id = int(parts[3])
                        cnt = int(parts[4])
                    except ValueError:
                        continue
                    counts[row_id] = cnt

            df[count_col] = df["ROW_ID"].map(lambda rid: counts.get(rid, 0)).astype(int)
            
            # If this is the ClinVar SV BED, also pull in the 5th BED column for phenotype
            if prefix == "clinvar_sv":
                clinvar_bed5_map = {}
                with open(raw_file, "r") as raw_fh2:
                    for line in raw_fh2:
                        parts = line.rstrip().split("\t")
                        if len(parts) < 9:
                            continue
                        try:
                            row_id = int(parts[3])
                        except ValueError:
                            continue
                        bed_col5 = parts[8]
                        if not bed_col5:
                            continue
                        if row_id not in clinvar_bed5_map:
                            clinvar_bed5_map[row_id] = set()
                        clinvar_bed5_map[row_id].add(bed_col5)

                def _join_bed5(rid: int) -> str:
                    vals = clinvar_bed5_map.get(rid)
                    if not vals:
                        return ""
                    return ";".join(sorted(vals))

                df["CLINVAR_PHENOTYPE"] = df["ROW_ID"].map(_join_bed5)

        except FileNotFoundError:
            print("[WARN] bedtools not found on PATH; skipping", count_col)
            df[count_col] = 0
            if prefix == "clinvar_sv":
                df["CLINVAR_PHENOTYPE"] = ""
        except subprocess.CalledProcessError as e:
            print(f"[WARN] bedtools intersect failed for {bed_path}: {e}")
            df[count_col] = 0
            if prefix == "clinvar_sv":
                df["CLINVAR_PHENOTYPE"] = ""

    # 1) Gene-level annotation: need both list and counts.
    if HG38_GENES_BED and os.path.exists(HG38_GENES_BED):
        genes_raw = os.path.join(out_dir, "hg38_genes_raw.tsv")
        try:
            with open(genes_raw, "w") as raw_fh:
                subprocess.run(
                    ["bedtools", "intersect", "-a", bed_input, "-b", HG38_GENES_BED, "-wa", "-wb"],
                    stdout=raw_fh,
                    check=True,
                )

            gene_hits = {}
            with open(genes_raw, "r") as gh:
                for line in gh:
                    parts = line.rstrip().split("\t")
                    if len(parts) < 8:
                        continue
                    try:
                        row_id = int(parts[3])
                    except ValueError:
                        continue
                    gene_name = parts[7]
                    if gene_name == ".":
                        continue
                    gene_hits.setdefault(row_id, set()).add(gene_name)

            df["GENE_LIST"] = df["ROW_ID"].map(
                lambda rid: ",".join(sorted(gene_hits.get(rid, set()))) if rid in gene_hits else ""
            )
            df["GENE_OVERLAP_COUNT"] = df["ROW_ID"].map(
                lambda rid: len(gene_hits.get(rid, set()))
            ).astype(int)

        except FileNotFoundError:
            print("[WARN] bedtools not found on PATH; skipping hg38 gene annotation.")
            df["GENE_LIST"] = ""
            df["GENE_OVERLAP_COUNT"] = 0
        except subprocess.CalledProcessError as e:
            print(f"[WARN] bedtools intersect failed for hg38 genes: {e}")
            df["GENE_LIST"] = ""
            df["GENE_OVERLAP_COUNT"] = 0
    else:
        print("[WARN] hg38 gene BED not found; skipping gene-level annotation.")
        df["GENE_LIST"] = ""
        df["GENE_OVERLAP_COUNT"] = 0

    # 2) Exon-level annotation (hg38 exons): counts only
    annotate_bed(HG38_EXONS_BED, "EXON_OVERLAP_COUNT", "hg38_exons")

    # 3) ClinGen HI / TS genes and recurrent CNV regions: counts only
    annotate_bed(CLINGEN_HI_GENES_BED, "CLINGEN_HI_OVERLAP_COUNT", "clingen_hi")
    annotate_bed(CLINGEN_TS_GENES_BED, "CLINGEN_TS_OVERLAP_COUNT", "clingen_ts")
    annotate_bed(CLINGEN_RECURRENT_CNV_BED, "CLINGEN_REC_CNV_OVERLAP_COUNT", "clingen_recurrent_cnv")
    
    # 4) ClinVar pathogenic / likely pathogenic SVs
    annotate_bed(CLINVAR_SV_BED, "CLINVAR_OVERLAP_COUNT", "clinvar_sv")

    return df


# 10. Pathogenicity-style flagging (Variants of Interest)

def flag_variants_of_interest(df: pd.DataFrame, child: str) -> pd.DataFrame:
    """
    Mark and keep variants that look clinically interesting.

    Uses simple tier labels (Tier0/1/2) based on ClinVar, ClinGen,
    size, exons and inheritance. Returns a filtered, sorted df.
    """

    df = df.copy()

    # Only consider variants where child is non-reference
    child_gt_col = f"{child}_GT"
    alt_gts = {"0/1", "1/0", "0|1", "1|0", "1/1", "1|1", "1"}
    if child_gt_col not in df.columns:
        raise ValueError(f"Child GT column not found: {child_gt_col}")
    df = df[df[child_gt_col].isin(alt_gts)].copy()


    def get_count(row, col):
        try:
            v = row.get(col, 0)
            if pd.isna(v):
                return 0
            return int(v)
        except Exception:
            return 0

    def classify(row):
        svtype = str(row.get("SVTYPE", "")).upper()
        size_cat = row.get("SIZE_PRIORITY", "unknown")
        inh = row.get("INHERITANCE", "unknown")

        exon_cnt   = get_count(row, "EXON_OVERLAP_COUNT")
        gene_cnt   = get_count(row, "GENE_OVERLAP_COUNT")
        hi_cnt     = get_count(row, "CLINGEN_HI_OVERLAP_COUNT")
        ts_cnt     = get_count(row, "CLINGEN_TS_OVERLAP_COUNT")
        rec_cnt    = get_count(row, "CLINGEN_REC_CNV_OVERLAP_COUNT")
        clinvar_cnt = get_count(row, "CLINVAR_OVERLAP_COUNT")
        
        # ----------------------------------------------------
        # TIER 0: ClinVar pathogenic / likely pathogenic SVs
        # ----------------------------------------------------
        if clinvar_cnt > 0:
            if inh == "de_novo":
                return "Tier0_clinvar_pathogenic_de_novo"
            else:
                return "Tier0_clinvar_pathogenic"


        # ---------------------------
        # TIER 1: strongest flags
        # ---------------------------

        # Recurrent CNV region 
        if rec_cnt > 0 and exon_cnt > 0 and svtype in {"DEL", "DUP"}:
            if inh == "de_novo":
                return "Tier1_recurrent_CNV_de_novo"
            else:
                return "Tier1_recurrent_CNV"

        # Deletions overlapping HI genes 
        if svtype == "DEL" and hi_cnt > 0 and exon_cnt > 0:
            if inh == "de_novo":
                return "Tier1_HI_del_de_novo"
            else:
                return "Tier1_HI_del"

        # Duplications overlapping TS genes
        if svtype == "DUP" and ts_cnt > 0 and exon_cnt > 0:
            if inh == "de_novo":
                return "Tier1_TS_dup_de_novo"
            else:
                return "Tier1_TS_dup"

        # ----------------------------------
        # TIER 2: strong but less specific
        # ----------------------------------

        # Large / moderate-size exonic events in the child
        if exon_cnt > 0 and size_cat in {"high", "moderate"}:
            if inh == "de_novo":
                return "Tier2_exonic_large_de_novo"
            else:
                return "Tier2_exonic_large"

        # -------------------------------------------
        # TIER 3: gene overlaps but weaker evidence
        # -------------------------------------------

        if gene_cnt > 0:
            return "Tier3_gene_overlap"

        # ---------------------------
        # TIER 4: everything else
        # ---------------------------
        return "Tier4_other"

    df["PATHOGENICITY_TIER"] = df.apply(classify, axis=1)

    # Keep only Tier 0, Tier 1, and Tier 2
    keep = (
        df["PATHOGENICITY_TIER"].str.startswith("Tier0") |
        df["PATHOGENICITY_TIER"].str.startswith("Tier1") |
        df["PATHOGENICITY_TIER"].str.startswith("Tier2")
    )

    df_voi = df[keep].copy()

    if "SVLEN_BP" in df_voi.columns:
        df_voi["SVLEN_BP"] = pd.to_numeric(df_voi["SVLEN_BP"], errors="coerce")
        df_voi.sort_values(
            by=["PATHOGENICITY_TIER", "SVLEN_BP"],
            ascending=[True, False],
            inplace=True
        )
    else:
        df_voi.sort_values(by=["PATHOGENICITY_TIER"], inplace=True)
        
    cols = list(df_voi.columns)
    if "PATHOGENICITY_TIER" in cols and "CLINVAR_PHENOTYPE" in cols:
        base_cols = [c for c in cols if c not in ("PATHOGENICITY_TIER", "CLINVAR_PHENOTYPE")]
        cols = base_cols + ["PATHOGENICITY_TIER", "CLINVAR_PHENOTYPE"]
        df_voi = df_voi[cols]
    elif "CLINVAR_PHENOTYPE" in cols:
        cols = [c for c in cols if c != "CLINVAR_PHENOTYPE"] + ["CLINVAR_PHENOTYPE"]
        df_voi = df_voi[cols]

    return df_voi



# 11. MAIN

def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)
    
    print("[INFO] Reading vcf...")

    # ---------------------------------------------
    # ASK USER WHICH 3 SAMPLES TO ANALYSE (by name)
    # ---------------------------------------------
    trio_samples = choose_trio_samples(args.vcf)

    print("\n[INFO] Selected samples:", ", ".join(trio_samples))
    print("[INFO] Starting initial AWK parsing...")

    # ------------------------------
    # RUN FIRST AWK PASS
    # ------------------------------
    parsed_tsv = os.path.join(args.out, "trio_table.tsv")

    print("[INFO] Running primary AWK parsing to extract trio genotype table...")

    with tempfile.NamedTemporaryFile("w", delete=False) as tf:
        tf.write(AWK_PARSE)
        awkfile = tf.name

    # pass chosen sample names into awk as variables sample1/2/3
    with open(parsed_tsv, "w") as o:
        subprocess.run(
            [
                "awk",
                "-v", f"sample1={trio_samples[0]}",
                "-v", f"sample2={trio_samples[1]}",
                "-v", f"sample3={trio_samples[2]}",
                "-f", awkfile,
                args.vcf,
            ],
            stdout=o,
            check=True,
        )


    os.unlink(awkfile)

    print("[INFO] Genotype table created:", parsed_tsv)

    # ------------------------------
    # LOAD PARSED TSV
    # ------------------------------
    print("[INFO] Loading parsed trio table...")

    df = pd.read_csv(parsed_tsv, sep="\t", dtype=str)

    # extract sample names from columns ending with _GT
    samples = [c.replace("_GT", "") for c in df.columns if c.endswith("_GT")]

    # ------------------------------
    # INFER TRIO ROLES + TRIO STATS
    # ------------------------------
    print("[INFO] Trio table loaded. Inferring trio roles...")

    roles, trio_stats = infer_trio(df, samples)
    print("[INFO] Determining who is child, father, mother...")

    child = roles["child"]
    p1 = roles["parent1"]
    p2 = roles["parent2"]

    min_vio = min(ts["violations"] for ts in trio_stats)

    print("\n===== TRIO ROLE SEARCH (Mendelian violations) =====")
    for ts in trio_stats:
        print(
            f"Assuming child={ts['child']}, "
            f"parent1={ts['parent1']}, parent2={ts['parent2']} "
            f"→ violations={ts['violations']}"
        )
    print(
        f"\nSelected trio = child={child}, parent1={p1}, parent2={p2}\n "
        f"because this assignment had the smallest number of Mendelian\n "
        f"violations ({min_vio}) among all permutations."
    )
    
    print("\n===== TRIO ROLES =====")
    print(f"Child    = {child}")
    print(f"Parent1  = {p1}")
    print(f"Parent2  = {p2}")

    # ---------------------------------------
    # INFER SEX FROM chrX  (child / parents)
    # ---------------------------------------
    print("\n[INFO] Performing sex inference from chrX heterozygosity...")

    sex_calls, x_het_rates = infer_sex_from_chrX(df, samples)

    child_sex = sex_calls.get(child, "unknown")
    p1_sex = sex_calls.get(p1, "unknown")
    p2_sex = sex_calls.get(p2, "unknown")

    father = None
    mother = None

    parent_sexes = {p1: p1_sex, p2: p2_sex}
    male_parents = [s for s, sx in parent_sexes.items() if sx == "M"]
    female_parents = [s for s, sx in parent_sexes.items() if sx == "F"]

    if male_parents:
        father = male_parents[0]
    if female_parents:
        mother = female_parents[0]

    print("\n===== SEX INFERENCE (chrX) =====")
    for s in samples:
        hr = x_het_rates.get(s, None)
        hr_str = "NA" if hr is None else f"{hr:.4f}"
        print(f"{s}: sex={sex_calls.get(s, 'unknown')}, X_het_rate={hr_str}")
    print(f"\nChild sex : {child} → {child_sex}")
    print(f"Father    : {father if father else 'undetermined'}")
    print(f"Mother    : {mother if mother else 'undetermined'}")

    # ------------------------------
    # FIND ORIGINAL COLUMN NUMBERS
    # ------------------------------
    header_samples = []
    with open(args.vcf) as f:
        for line in f:
            if line.startswith("#CHROM"):
                header_samples = line.strip().split("\t")
                break

    sample_to_col = {}
    for idx, name in enumerate(header_samples):
        if idx >= 9:
            sample_to_col[name] = idx + 1 

    child_col = sample_to_col[child]
    p1_col = sample_to_col[father]
    p2_col = sample_to_col[mother]

    print("\n===== VCF COLUMN NUMBERS =====")
    print(f"Child col   = {child_col}")
    print(f"Father col = {p1_col}")
    print(f"Mother col = {p2_col}")
    
    # ------------------------------
    # RUN AWK DE-NOVO DETECTOR
    # ------------------------------
    print("\n[INFO] Running de novo variant discovery ...")

    denovo_file_txt = os.path.join(args.out, "de_novo_calls.txt")
    denovo_file_tsv = os.path.join(args.out, "de_novo_calls.tsv")
    run_awk_denovo(args.vcf, denovo_file_txt, denovo_file_tsv, p1_col, p2_col, child_col)
    print("[INFO] De novo variant detection completed.")
    
    # ------------------------------------------------
    # ANNOTATE SV LENGTH, SIZE PRIORITY & INHERITANCE
    # ------------------------------------------------
    print("[INFO] Getting SV Type Counts ...")
    print("[INFO] Annotating variants with ClinVar & ClinGen datasets...")

    df = add_sv_length_and_size_category(df)
    df = annotate_inheritance(df, child=child, p1=p1, p2=p2)

    # Genome annotation via bedtools (hg38 genes/exons + ClinGen dosage/CNVs) on ALL SVs
    df = bedtools_annotate_all_variants(df, args.out)

    # child non-ref subset (these are the variants actually present in child)
    child_gt_col = f"{child}_GT"
    alt_gts = {"0/1", "1/0", "0|1", "1|0", "1/1", "1|1", "1"}
    df_child_alt = df[df[child_gt_col].isin(alt_gts)].copy()
    
    # ------------------------------
    # COUNT ALL SVTYPES
    # ------------------------------
    sv_counts = df["SVTYPE"].value_counts().reset_index()
    sv_counts.columns = ["SVTYPE", "COUNT"]
    sv_counts.to_csv(os.path.join(args.out, "svtype_counts.tsv"), sep="\t", index=False)
    
    print("\n===== SVTYPE COUNTS =====")
    print(str(sv_counts))
    
    # ------------------------------------------
    # Pathogenicity-style variants of interest
    # ------------------------------------------
    print("\n[INFO] Annotation completed. Filtering variants of interest")

    df_voi = flag_variants_of_interest(df, child=child)
    voi_path = os.path.join(args.out, "clinically_prioritised_SVs.tsv")
    df_voi.to_csv(voi_path, sep="\t", index=False)
    print(f"\nClinically prioritised SVs (Tier 0/1/2) written to: {voi_path}")

    if not df_voi.empty:
        voi_tier_counts = df_voi["PATHOGENICITY_TIER"].value_counts().reset_index()
        voi_tier_counts.columns = ["PATHOGENICITY_TIER", "COUNT"]
    else:
        voi_tier_counts = pd.DataFrame(columns=["PATHOGENICITY_TIER", "COUNT"])
    
    # ------------------------------
    # Overall Counts from Bedtools
    # ------------------------------
    total_sv = len(df)

    gene_ov_sv = int((df["GENE_OVERLAP_COUNT"] > 0).sum()) if "GENE_OVERLAP_COUNT" in df.columns else None
    exon_ov_sv = int((df["EXON_OVERLAP_COUNT"] > 0).sum()) if "EXON_OVERLAP_COUNT" in df.columns else None
    hi_ov_sv = int((df["CLINGEN_HI_OVERLAP_COUNT"] > 0).sum()) if "CLINGEN_HI_OVERLAP_COUNT" in df.columns else None
    ts_ov_sv = int((df["CLINGEN_TS_OVERLAP_COUNT"] > 0).sum()) if "CLINGEN_TS_OVERLAP_COUNT" in df.columns else None
    rec_cnv_ov_sv = int((df["CLINGEN_REC_CNV_OVERLAP_COUNT"] > 0).sum()) if "CLINGEN_REC_CNV_OVERLAP_COUNT" in df.columns else None
    clinvar_ov_sv = int((df["CLINVAR_OVERLAP_COUNT"] > 0).sum()) if "CLINVAR_OVERLAP_COUNT" in df.columns else None


    # ------------------------------
    # REPORT
    # ------------------------------
    rep = os.path.join(args.out, "report.txt")
    with open(rep, "w") as f:

        f.write("===== TRIO ROLE SEARCH (Mendelian violations) =====\n")
        for ts in trio_stats:
            f.write(
                f"Assuming child={ts['child']}, "
                f"parent1={ts['parent1']}, parent2={ts['parent2']} "
                f"→ violations={ts['violations']}\n"
            )
        f.write("\n")
        f.write(
            f"Selected trio = child={child}, parent1={p1}, parent2={p2} because this \n"
            f"assignment had the smallest number of Mendelian violations \n"
            f"({min_vio}) among all permutations.\n\n"
        )
        
        f.write("===== TRIO ROLES =====\n")
        f.write(f"Child    : {child}\n")
        f.write(f"Parent1  : {p1}\n")
        f.write(f"Parent2  : {p2}\n\n")

        f.write("===== COLUMN NUMBERS =====\n")
        f.write(f"Child col   : {child_col}\n")
        f.write(f"Parent1 col : {p1_col}\n")
        f.write(f"Parent2 col : {p2_col}\n\n")


        f.write("===== SVTYPE COUNTS =====\n")
        f.write(str(sv_counts))
        f.write("\n\n")
        f.write("===== De Novo Calls =====\n")
        f.write(f"Written to: {denovo_file_txt} and {denovo_file_tsv}\n\n")
        
        f.write("=== De Novo Variants ===\n")
        def extract_gt(sample_field: str) -> str:
            """Return GT from a VCF sample field like '0/1:35:.'."""
            sample_field = sample_field.strip()
            if sample_field in ("", ".", "./."):
                return sample_field or "./."
            return sample_field.split(":", 1)[0]

        with open(denovo_file_txt, "r") as tsv_fh:
            first_data_line = True
            for line in tsv_fh:
                line = line.strip()
                if not line:
                    continue

                if first_data_line:
                    if line.startswith("CHROM\tSTART_POS\tEND_POS\tID\tREF\tALT\tSV_TYPE"):
                        first_data_line = False
                        continue
                    first_data_line = False

                if line.startswith("#"):
                    continue

                cols = line.split("\t")
                if len(cols) < 12:
                    continue

                chrom  = cols[0]
                pos    = cols[1]
                _id    = cols[2]
                ref    = cols[3]
                alt    = cols[4]
                qual   = cols[5]
                flt    = cols[6]
                info   = cols[7]
                fmt    = cols[8]
                p1_sf  = cols[9]
                p2_sf  = cols[10]
                ch_sf  = cols[11]

                p1_gt    = extract_gt(p1_sf)
                p2_gt    = extract_gt(p2_sf)
                child_gt = extract_gt(ch_sf)

                svtype = "NA"
                for item in info.split(";"):
                    if item.startswith("SVTYPE="):
                        svtype = item.split("=", 1)[1]
                        break

                f.write(
                    f"{chrom}:{pos}  "
                    f"REF={ref}  ALT={alt}  QUAL={qual}  FILTER={flt}  SVTYPE={svtype}  "
                    f"{child}_GT={child_gt}, "
                    f"{p1}_GT={p1_gt}, "
                    f"{p2}_GT={p2_gt}\n"
                )
        
        f.write("\n\n")
        
        f.write("===== ANNOTATION using bedtools + hg38 / Clinvar / ClinGen =====\n")
        f.write(f"Total SVs analysed                       : {total_sv}\n")
        f.write(f"hg38 gene annotation BED                 : {HG38_GENES_BED}\n")
        f.write(f"hg38 exon annotation BED                 : {HG38_EXONS_BED}\n")
        f.write(f"ClinGen HI genes BED (GRCh38)            : {CLINGEN_HI_GENES_BED}\n")
        f.write(f"ClinGen TS genes BED (GRCh38)            : {CLINGEN_TS_GENES_BED}\n")
        f.write(f"ClinGen recurrent CNV BED (hg38)         : {CLINGEN_RECURRENT_CNV_BED}\n")
        f.write(f"ClinVar pathogenic SV BED (GRCh38)       : {CLINVAR_SV_BED}\n")
        f.write("\n")
        if gene_ov_sv is not None:
            f.write(f"SVs overlapping ≥1 hg38 gene             : {gene_ov_sv} / {total_sv}\n")
        if exon_ov_sv is not None:
            f.write(f"SVs overlapping ≥1 hg38 exon             : {exon_ov_sv} / {total_sv}\n")
        if hi_ov_sv is not None:
            f.write(f"SVs overlapping ≥1 ClinGen HI gene       : {hi_ov_sv} / {total_sv}\n")
        if ts_ov_sv is not None:
            f.write(f"SVs overlapping ≥1 ClinGen TS gene       : {ts_ov_sv} / {total_sv}\n")
        if rec_cnv_ov_sv is not None:
            f.write(f"SVs overlapping ≥1 recurrent CNV region  : {rec_cnv_ov_sv} / {total_sv}\n")
        if clinvar_ov_sv is not None:
            f.write(f"SVs overlapping ≥1 ClinVar pathogenic SV : {clinvar_ov_sv} / {total_sv}\n")
        f.write("\n\n")
        
        f.write("===== CLINICALLY PRIORITISED SVs based on PATHOGENICITY TIERS =====\n")
        f.write(f"Output file: {os.path.basename(voi_path)}\n")
        if df_voi.empty:
            f.write("No Tier1 or Tier2 variants of interest were identified.\n\n")
        else:
            f.write("Tier counts (Tier 0/1/2 only):\n")
            for _, row in voi_tier_counts.iterrows():
                f.write(f"  {row['PATHOGENICITY_TIER']:30s} : {row['COUNT']}\n")
            f.write("\n")
            f.write("Top candidate variants (first 10):\n")
            cols_basic = ["CHROM", "POS", "END_INT", "SVTYPE", "SVLEN_BP",
                          "GENE_LIST", "INHERITANCE", "PATHOGENICITY_TIER", "CLINVAR_PHENOTYPE"]
            cols_basic = [c for c in cols_basic if c in df_voi.columns]
            head10 = df_voi[cols_basic].head(10)

            f.write(head10.to_string(index=False))
            f.write("\n\n")

    print(f"\nReport saved in: {rep}")
    
    
    # Inheritance category counts
    inheritance_counts = df["INHERITANCE"].value_counts(dropna=False).reset_index()
    inheritance_counts.columns = ["INHERITANCE", "COUNT"]
    
    inh_tsv = os.path.join(args.out, "inheritance_counts.tsv")
    inheritance_counts.to_csv(inh_tsv, sep="\t", index=False)
    

    # ----------------------------------------------
    # Generate plots using custom plotting package
    # ----------------------------------------------
    plot_dir = os.path.join(args.out, "plots")
    os.makedirs(plot_dir, exist_ok=True)
    make_plots(plot_dir)
    
    # ------------------------------
    # FULL SUMMARY
    # ------------------------------
    summary_path = os.path.join(args.out, "summary_full.txt")

    
    with open(summary_path, "w") as f2:
        f2.write("VCF Trio Structural Variant Analysis – Full Summary\n")
        f2.write("==================================================\n\n")

        # 1. Trio composition, Mendelian search and sex
        f2.write("1. Trio Composition, Mendelian Search and Sex Inference\n")

        f2.write("   1.1 Trio role search (Mendelian violations across all permutations)\n")
        for ts in trio_stats:
            f2.write(
                f"       - child={ts['child']}, "
                f"parent1={ts['parent1']}, parent2={ts['parent2']} "
                f"→ violations={ts['violations']}\n"
            )
        f2.write(
            f"\n       Selected trio = child={child}, parent1={p1}, parent2={p2} \n"
            f"        because this assignment had the smallest number of Mendelian\n "
            f"        violations ({min_vio}) among all permutations.\n\n"
        )

        f2.write("   1.2 Trio members and sex calls\n")
        f2.write(f"       - Child   : {child} (sex={child_sex})\n")
        f2.write(f"       - Parent1 : {p1} (sex={p1_sex})\n")
        f2.write(f"       - Parent2 : {p2} (sex={p2_sex})\n")
        f2.write(f"       - Father  : {father or 'undetermined'}\n")
        f2.write(f"       - Mother  : {mother or 'undetermined'}\n\n")
        f2.write("       chrX heterozygosity rates used for sex calling:\n")
        for s in samples:
            hr = x_het_rates.get(s, None)
            hr_str = "NA" if hr is None else f"{hr:.4f}"
            f2.write(
                f"         * {s:10s} : "
                f"sex={sex_calls.get(s, 'unknown')}, X_het_rate={hr_str}\n"
            )
        f2.write("\n")

        # 2. Global SV statistics
        f2.write("2. Global Structural Variant Counts (all samples, all SVs)\n")
        f2.write(sv_counts.to_string(index=False))
        f2.write("\n\n")

        # 3. Inheritance categories
        f2.write("3. Child-Centric Inheritance Classification\n")
        f2.write("   Each SV is labelled based on child + parents' genotypes.\n")
        f2.write("   Categories: de_novo, inherited_p1, inherited_p2, inherited_both,\n")
        f2.write("   reference_or_absent_child, unknown.\n\n")
        f2.write(inheritance_counts.to_string(index=False))
        f2.write("\n\n")
        
        # 4. De Novo Variants
        f2.write("4. De Novo Variants\n")
        def extract_gt(sample_field: str) -> str:
            """Return GT from a VCF sample field like '0/1:35:.'."""
            sample_field = sample_field.strip()
            if sample_field in ("", ".", "./."):
                return sample_field or "./."
            return sample_field.split(":", 1)[0]

        with open(denovo_file_txt, "r") as tsv_fh:
            first_data_line = True
            for line in tsv_fh:
                line = line.strip()
                if not line:
                    continue

                if first_data_line:
                    if line.startswith("CHROM\tSTART_POS\tEND_POS\tID\tREF\tALT\tSV_TYPE"):
                        first_data_line = False
                        continue
                    first_data_line = False

                if line.startswith("#"):
                    continue

                cols = line.split("\t")
                if len(cols) < 12:
                    continue

                chrom  = cols[0]
                pos    = cols[1]
                _id    = cols[2]
                ref    = cols[3]
                alt    = cols[4]
                qual   = cols[5]
                flt    = cols[6]
                info   = cols[7]
                fmt    = cols[8]
                p1_sf  = cols[9]
                p2_sf  = cols[10]
                ch_sf  = cols[11]

                # get GTs from sample fields
                p1_gt    = extract_gt(p1_sf)
                p2_gt    = extract_gt(p2_sf)
                child_gt = extract_gt(ch_sf)

                # parse SVTYPE from INFO if present
                svtype = "NA"
                for item in info.split(";"):
                    if item.startswith("SVTYPE="):
                        svtype = item.split("=", 1)[1]
                        break

                f2.write(
                    f"{chrom}:{pos}  "
                    f"REF={ref}  ALT={alt}  QUAL={qual}  FILTER={flt}  SVTYPE={svtype}  "
                    f"{child}_GT={child_gt}, "
                    f"{p1}_GT={p1_gt}, "
                    f"{p2}_GT={p2_gt}\n"
                )
        
        f2.write("\n\n")
        
        # 5. Genome annotation (hg38 genes + ClinVar/ClinGen)
        f2.write("5. Genome Annotation using hg38 genes, ClinVar/ClinGen\n")
        f2.write("   Annotation performed with bedtools intersect using:\n")
        f2.write(f"     - hg38 gene annotation BED                 : {HG38_GENES_BED}\n")
        f2.write(f"     - hg38 exon annotation BED                 : {HG38_EXONS_BED}\n")
        f2.write(f"     - ClinGen HI genes BED (GRCh38)           : {CLINGEN_HI_GENES_BED}\n")
        f2.write(f"     - ClinGen TS genes BED (GRCh38)           : {CLINGEN_TS_GENES_BED}\n")
        f2.write(f"     - ClinGen recurrent CNV BED (hg38)        : {CLINGEN_RECURRENT_CNV_BED}\n")
        f2.write(f"     - ClinVar pathogenic SV BED (GRCh38)      : {CLINVAR_SV_BED}\n")
        f2.write("   Summary of SVs with ≥1 overlap:\n")
        if gene_ov_sv is not None:
            f2.write(f"     * SVs overlapping ≥1 hg38 gene            : {gene_ov_sv} / {total_sv}\n")
        if exon_ov_sv is not None:
            f2.write(f"     * SVs overlapping ≥1 hg38 exon            : {exon_ov_sv} / {total_sv}\n")
        if hi_ov_sv is not None:
            f2.write(f"     * SVs overlapping ≥1 ClinGen HI gene      : {hi_ov_sv} / {total_sv}\n")
        if ts_ov_sv is not None:
            f2.write(f"     * SVs overlapping ≥1 ClinGen TS gene      : {ts_ov_sv} / {total_sv}\n")
        if rec_cnv_ov_sv is not None:
            f2.write(f"     * SVs overlapping ≥1 recurrent CNV region : {rec_cnv_ov_sv} / {total_sv}\n")
        if clinvar_ov_sv is not None:
            f2.write(f"     * SVs overlapping ≥1 ClinVar pathogenic SV : {clinvar_ov_sv} / {total_sv}\n")
        f2.write("\n")

        
        # 6. Clinically prioritised SVs
        f2.write("6. Clinically Prioritised SVs (VOI – Variants of Interest)\n")
        f2.write("   Using existing annotations, each SV present in the child is assigned\n")
        f2.write("   a PATHOGENICITY_TIER based on:\n")
        f2.write("     - SVTYPE (DEL, DUP, INS, etc.)\n")
        f2.write("     - Size category (SIZE_PRIORITY: high / moderate / low)\n")
        f2.write("     - Inheritance (INHERITANCE: de_novo, inherited_p1, etc.)\n")
        f2.write("     - Exon overlap (EXON_OVERLAP_COUNT)\n")
        f2.write("     - Overlap with ClinGen HI / TS genes and recurrent CNV regions.\n")
        f2.write("\n")
        f2.write("   Tier definitions (heuristic):\n")
        f2.write("     - Tier0_clinvar_pathogenic* : SVs overlapping ClinVar Pathogenic/\n")
        f2.write("                                   Likely pathogenic CNVs (GRCh38);\n")
        f2.write("                                   de novo subclass flagged.\n")
        f2.write("     - Tier1_recurrent_CNV* : SV overlapping a ClinGen recurrent CNV\n")
        f2.write("                              region (DEL/DUP, exonic), often known\n")
        f2.write("                              pathogenic microdeletion/duplication\n")
        f2.write("                              syndromes; de novo subclass flagged.\n")
        f2.write("     - Tier1_HI_del*        : Deletions overlapping ClinGen HI genes\n")
        f2.write("                              in exons; de novo subclass flagged.\n")
        f2.write("     - Tier1_TS_dup*        : Duplications overlapping ClinGen TS genes\n")
        f2.write("                              in exons; de novo subclass flagged.\n")
        f2.write("     - Tier2_exonic_large*  : Large/moderate exonic SVs without direct\n")
        f2.write("                              HI/TS evidence but potentially impactful.\n")
        f2.write("     - Tier3_gene_overlap   : SVs overlapping genes but not exons/\n")
        f2.write("                              dosage-sensitive categories.\n")
        f2.write("     - Tier4_other          : Remaining SVs.\n\n")

        if df_voi.empty:
            f2.write("   No Tier1 or Tier2 variants of interest were identified.\n\n")
        else:
            f2.write("   Counts of Tier 0/1/2 variants:\n")
            for _, row in voi_tier_counts.iterrows():
                f2.write(f"     * {row['PATHOGENICITY_TIER']:30s} : {row['COUNT']}\n")
            f2.write("\n")
            f2.write("   A detailed table of Tier 0/1/2 variants is written to:\n")
            f2.write(f"     - {os.path.basename(voi_path)}\n\n")

            
        # 7. Output files guide
        f2.write("7. Output Files\n")
        f2.write("   - trio_table.tsv              : Parsed AWK table with GTs and SV annotations.\n")
        f2.write("   - svtype_counts.tsv           : Global count of SVTYPEs.\n")
        f2.write("   - de_novo_calls.txt / .tsv    : De novo candidate SVs (0/0 parents, child non-ref).\n")
        f2.write("   - clinically_prioritised_SVs.tsv : Child's SVs passing size/region filters.\n")
        f2.write("   - sv_for_bedtools.bed         : BED intervals used for all bedtools annotations.\n")
        f2.write("   - bedtools_hg38_genes_raw.tsv : Raw SV × gene overlaps.\n")
        f2.write("   - hg38_exons_raw.tsv / _counts.tsv : Raw + counts for exon overlaps.\n")
        f2.write("   - clingen_hi_raw.tsv / _counts.tsv : Raw + counts for ClinGen HI genes.\n")
        f2.write("   - clingen_ts_raw.tsv / _counts.tsv : Raw + counts for ClinGen TS genes.\n")
        f2.write("   - clingen_recurrent_cnv_raw.tsv / _counts.tsv : Raw + counts for recurrent CNVs.\n")
        f2.write("   - report.txt                  : Compact summary (key results).\n")
        f2.write("   - summary_full.txt            : Full overview of findings.\n")
        
        f2.write("\n")
        f2.write("8. Visualisation Outputs (Figures)\n")
        f2.write("     - fig1_svtype_counts.png\n")
        f2.write("         Global distribution of structural variant types (DEL, INS, DUP).\n")
        f2.write("     - fig2_inheritance_counts.png\n")
        f2.write("         Inheritance pattern of SVs in the child, showing de novo and\n")
        f2.write("         inherited categories (inherited_p1, inherited_p2, inherited_both,\n")
        f2.write("         reference_or_absent_child, unknown).\n")
        f2.write("     - fig3_annotation_overlap.png\n")
        f2.write("         Number of SVs with ≥1 overlap with hg38 exons, ClinGen\n")
        f2.write("         haploinsufficiency (HI) genes, triplosensitivity (TS) genes,\n")
        f2.write("         recurrent CNV regions, and ClinVar pathogenic structural variants.\n")
        f2.write("     - fig4_pathogenicity_tiers.png\n")
        f2.write("         Distribution of clinically prioritised SVs by PATHOGENICITY_TIER\n")
        f2.write("         (e.g. Tier0_clinvar_pathogenic, Tier1_HI_del, Tier1_recurrent_CNV,\n")
        f2.write("         Tier2_exonic_large, Tier2_exonic_large_de_novo).\n")
        f2.write("     - fig5_per_chrom_svtype.png\n")
        f2.write("         Overview of number of SV type per chromosome.\n")
        f2.write("\n")


    print(f"Full summary saved in   : {summary_path}")


    print("\n[INFO] Analysis completed successfully. Results available in:", args.out)
    print("[INFO] Thank you for using the trio analyzer!\n")


if __name__ == "__main__":
    main()
