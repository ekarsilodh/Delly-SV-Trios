# Usage Notes for VCF Trio Analyzer

This document contains extended notes, tips, and practical guidance for using
and customising the `vcf_analyzer.py` pipeline. For general instructions, see
the main `README.md`. This file is intended for users who want additional
technical details or want to extend the tool.

---

## 1. Overview of the Pipeline

The pipeline performs the following major steps:

1. **VCF Header Parsing**
   - Reads sample names from the `#CHROM` header line.
   - Lets the user manually select 3 samples to form a trio.

2. **Trio Role Inference**
   - Tests all permutations of the 3 samples.
   - Calculates Mendelian violations for each combination.
   - Selects: **child, mother, father** configuration with fewest violations.

3. **Sex Inference (chrX Heterozygosity)**
   - For each sample:
     - Computes ratio of heterozygous genotypes on chrX.
   - Classifies:
     - Low heterozygosity → **Male**
     - High heterozygosity → **Female**

4. **AWK-Based Parsing**
   - A fast AWK command processes the raw structural variant VCF.
   - Extracts: chrom, pos, end, SVTYPE, gene fields, and trio genotypes.

5. **De Novo Detection**
   - Simple rule:
     - Both parents `0/0`
     - Child ≠ `0/0`
   - Writes results to:
     - `de_novo_calls.txt`
     - `de_novo_calls.tsv`

6. **Annotation (via bedtools)**
   - Intersect SV intervals with:
     - hg38 genes
     - hg38 exons
     - ClinGen HI/TS genes
     - ClinGen recurrent CNVs
     - ClinVar pathogenic SVs

7. **Tiering / Prioritisation**
   - Uses size, region, and annotation features to classify variants into:
     - Tier 0 (highest importance)
     - Tier 1
     - Tier 2

8. **Report Generation**
   - Creates `report.txt` and `summary_full.txt`.

9. **Optional Plotting**
   - If `make_plots(out_dir)` is implemented in `src/plot.py`,
     visualisations are saved inside the output directory.

---

## 2. Running the Pipeline: Basic Example

```bash
python src/vcf_analyzer.py \
  --vcf path/to/multi_sample.vcf \
  --out out/

## 3. Expected Output Files

Inside your output directory:

| File | Description |
|------|-------------|
| `trio_table.tsv` | Parsed SV table with trio genotype columns |
| `svtype_counts.tsv` | Count per SVTYPE |
| `de_novo_calls.txt` / `de_novo_calls.tsv` | De novo candidates |
| `clinically_prioritised_SVs.tsv` | Tiered child variants |
| `report.txt` | Short summary |
| `summary_full.txt` | Detailed multi-section report |
| `*.png` | Optional plots |

---

## 4. Customising the Pipeline

### 4.1 Changing BED File Paths

At the top of `src/vcf_analyzer.py`, update:

```python
HG38_GENES_BED = "hg38_genes.bed"
HG38_EXONS_BED = "hg38_exons.bed"
CLINGEN_HI_GENES_BED = "ClinGen_haploinsufficiency_gene_GRCh38.bed"
CLINGEN_TS_GENES_BED = "ClinGen_triplosensitivity_gene_GRCh38.bed"
CLINGEN_RECURRENT_CNV_BED = "ClinGen recurrent CNV-hg38.bed"
CLINVAR_SV_BED = "clinvar_SV_clean.bed"
```

Change them to absolute or relative paths, e.g.:

```python
HG38_GENES_BED = "databases/hg38_genes.bed"
```

---

### 4.2 Modifying Sex Inference Thresholds

Inside `infer_sex_from_chrX()` the defaults are:

```python
het_threshold = 0.1
min_sites = 5
```

Increase `min_sites` for very small SV VCFs to avoid unstable estimates.

---

### 4.3 Adjusting De Novo Logic

The default de novo detection rule is:

- `Parent1 = 0/0`
- `Parent2 = 0/0`
- `Child != 0/0`

You may extend this to include genotype quality (GQ), depth (DP), or filtering (`FILTER` field).

---

## 5. Troubleshooting

### ❗ `bedtools: command not found`

Install bedtools:

```bash
brew install bedtools         # macOS (Homebrew)
sudo apt install bedtools     # Ubuntu/Debian
```

---

### ❗ Output files missing

Check:

- BED files exist
- Paths to BED files are correct
- The VCF contains `SVTYPE` annotations

---

### ❗ Genotypes appear as missing (".")

This occurs if:

- Your VCF uses a different FORMAT order
- FORMAT fields like GT are missing for some samples

---

## 6. Recommended Practice

- Ensure your VCF and BED files use **consistent genome builds** (GRCh38 only).
- Use BGZIP + TABIX-indexed VCFs for best performance.
- Keep BED files sorted (`sort -k1,1 -k2,2n`).

---

✔ Ready to run the pipeline.