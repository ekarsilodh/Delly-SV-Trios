# Delly VCF Trio Analyzer

[![DOI](https://zenodo.org/badge/1106269273.svg)](https://doi.org/10.5281/zenodo.17781904)

A lightweight, research-oriented pipeline for **trio-based structural variant (SV) analysis** using a combination of Python, AWK, and BEDTools.  
The tool is designed to work on **multi-sample Delly VCF (GRCh38/hg38)** files containing structural variants.

<p align="center">
  <img src="banner_dark.png" width="900">
</p>

This pipeline:

- Allows users to manually choose **three samples** to form a trio
- Automatically infers **child / mother / father** roles using Mendelian violations
- Computes **X chromosome heterozygosity** for sex inference
- Performs fast **AWK genotype parsing**
- Detects **de novo structural variants**
- Annotates SVs using **bedtools intersect** with multiple BED databases
- Prioritizes variants into **Tier 0/1/2**
- Generates summary TSV files, text reports, and optional plots

> ⚠️ **For research and teaching only. Not a clinical diagnostic pipeline.**  

---

# 📦 Repository Structure

```
DELLY-SV-Trios/
├── README.md
├── LICENSE
├── requirements.txt
├── .gitignore
│
├── src/
│   ├── vcf_analyzer.py      # main analysis script
│   └── plot.py              # optional plotting functions
│
├── databases/
│   └── README.md            # explains required BED annotation files
│
├── examples/
│   └── example_run.sh       # advanced example run script
│
├── docs/
│   └── usage_notes.md       # extended notes, customisation tips
│
└── out/                     # output directory (ignored by git)
```

---

# 🛠️ Software Dependencies

The following software must be installed:

| Dependency | Purpose |
|-----------|---------|
| **Python ≥ 3.8** | Runs the main script |
| **awk** | Fast VCF parsing |
| **bedtools** | Annotation using genomic intervals |
| **git** | Version control (optional) |

### Check installations

```bash
python3 --version
awk --version
bedtools --version
```

---

# 📦 Python Dependencies

Listed in `requirements.txt`:

```
pandas
matplotlib
```

Install via:

```bash
pip install -r requirements.txt
```

---

# 📥 Required Input Files

## 1️⃣ Multi-sample VCF

Must be:

- GRCh38/hg38
- Structural variant VCF
- Containing three or more samples

## 2️⃣ BED Files (Mandatory)

- hg38_genes.bed  
- hg38_exons.bed  
- ClinGen_haploinsufficiency_gene_GRCh38.bed  
- ClinGen_triplosensitivity_gene_GRCh38.bed  
- ClinGen recurrent CNV-hg38.bed  
- clinvar_SV_clean.bed  

Full details: `databases/README.md`

---

# 🚀 Usage

```bash
python src/vcf_analyzer.py --vcf path/to/trio.vcf --out out/
```

---

# 🧪 Example Script

See:

```
examples/example_run.sh
```

Run:

```bash
chmod +x examples/example_run.sh
./examples/example_run.sh
```

---

# 📊 Output Files

- trio_table.tsv  
- svtype_counts.tsv  
- de_novo_calls.tsv  
- clinically_prioritised_SVs.tsv  
- report.txt  
- summary_full.txt  
- PNG plots (optional)

---

# 🔧 Customisation

Edit BED paths in `vcf_analyzer.py`:

```python
HG38_GENES_BED = "databases/hg38_genes.bed"
```

---

# 📚 Documentation

Additional notes:

```
docs/usage_notes.md
```

---

# 📝 License

MIT License (see LICENSE)

---

## If this project was helpful and you have used this in your work, please cite the following:
E. Lodh, S. Panda, D. Kadrekar “Delly-VCF-Trio-Analyzer: A Complete Pipeline for Trio-Based Structural Variant Interpretation”. Zenodo, Dec. 01, 2025. doi: [10.5281/zenodo.17781905](https://doi.org/10.5281/zenodo.17781905).

