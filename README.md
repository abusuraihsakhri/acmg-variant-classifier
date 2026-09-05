# ACMG Variant Classifier & Functional Genomics Annotator

A Python clinical genetics and bioinformatics library and CLI tool. Evaluates germline sequence variants using the 2015 ACMG/AMP 28-evidence-code guideline combination rules, coupled with the Tavtigian et al. (2018/2020) Bayesian point-scoring and posterior probability model. Also provides automated population allele frequency checking (gnomAD / 1000 Genomes thresholds), position-weight splice site impact prediction, CPIC pharmacogenomic variant annotation, and protein functional domain mapping.

Requires Python standard library only (zero external runtime dependencies).

---

## Features

- **ACMG/AMP 2015 Categorical Classification:** Implements the complete Table 5 rule combinations (PVS1, PS1-4, PM1-6, PP1-5, BA1, BS1-4, BP1-7) categorizing variants into Pathogenic, Likely Pathogenic, Uncertain Significance (VUS), Likely Benign, or Benign.
- **Tavtigian 2020 Bayesian Scoring & Posterior Probability:** Maps evidence tiers onto a continuous point scale (Supporting=1, Moderate=2, Strong=4, Very Strong=8) to compute cumulative points and calculate Bayesian posterior probabilities $P(\text{Pathogenicity} \mid \text{Evidence})$.
- **Automated Allele Frequency Evaluation:** Cross-checks population allele frequency against stand-alone benign ($BA1 > 0.05$) and absent/rare ($PM2 \le 0.0001$) thresholds, detecting contradictions.
- **In-Silico Splice Site Disruption:** Models canonical donor/acceptor site loss, cryptic activation, and exonic splice alterations.
- **PharmacoGenomics (CPIC):** Annotates Level A/B clinical pharmacogenomics guidelines across actionable genes (CYP2D6, CYP2C19, CYP2C9, DPYD, TPMT, NUDT15, SLCO1B1, VKORC1, HLA-B, IFNL3).
- **Functional Domain & Hotspot Mapping:** Cross-references variant amino acid positions against known protein domains and mutational hotspots (BRCA1, TP53, CFTR, BRAF, KRAS, PIK3CA).
- **Multiple Input Modes:** Interactive CLI wizard, single-variant flags, batch CSV/TSV table evaluation, and annotated VCF processing.

---

## Installation & Requirements

- Python 3.10+ (tested on 3.10, 3.11, 3.12)
- Zero external runtime dependencies. `pytest` is optional for running the test suite.

```bash
git clone https://github.com/abusuraihsakhri/acmg-variant-classifier.git
cd acmg-variant-classifier
```

## Security & Input Validation

The CLI includes built-in input validation:

- **Path validation**: Input files are verified to exist and be regular files before processing. Output paths are validated to prevent accidental overwrites of directories.
- **Probability range checking**: Allele frequency (`--af`), BA1 threshold (`--ba1-threshold`), and PM2 threshold (`--pm2-threshold`) values are validated to be within the [0, 1] range.
- **Graceful error handling**: File I/O errors and invalid inputs produce clear error messages with non-zero exit codes.

---

## CLI Usage

### 1. Single Variant Evaluation
Evaluate evidence codes directly:
```bash
python cli.py --evidence PVS1 PS3 PM2 --variant-id "BRCA1:c.5266dupC"
```
Output as JSON:
```bash
python cli.py --evidence PVS1 PS3 PM2 --json
```

Evaluate with population allele frequency:
```bash
python cli.py --evidence PS1 PS3 PM1 PP3 --af 0.000005 --variant-id "TP53:c.524G>A" --json
```

### 2. Batch Processing (CSV / TSV)
Run evaluation over a tabular file with `variant_id`, `evidence`, and optional `af`:
```bash
python cli.py -i sample_variants.csv --format json
```

### 3. VCF Processing
Parse VCF files containing `ACMG=` and optional `AF=` tags in INFO fields:
```bash
python cli.py --vcf variants.vcf --output results.txt
```

### 4. Splice Impact Prediction
Predict splice donor/acceptor disruption:
```bash
python cli.py --splice --pos 1 --ref G --alt A --distance 1 --json
```

### 5. Pharmacogenomics (CPIC) Annotation
Annotate actionable drug-gene guidelines:
```bash
python cli.py --pgx --gene CYP2C19 --variant-id "c.681G>A" --json
```

### 6. Protein Domain & Hotspot Mapping
Map variant position to functional domains:
```bash
python cli.py --domain --gene TP53 --pos 175 --json
```

### 7. Interactive Evaluation Wizard
Step-by-step guidance in the terminal:
```bash
python cli.py --interactive
```

---

## Python API Quickstart

```python
from acmg_classifier import classify_variant, report_to_text, report_to_dict
from variant_analysis import SpliceImpactPredictor, PharmacoGenomicsAnnotator, FunctionalDomainMapper

# 1. Classify a variant
report = classify_variant(
    codes=["PVS1", "PS3", "PM2"],
    variant_id="BRCA1:c.5266dupC",
    population_af=0.00002
)
print("Classification:", report.final_classification)  # Pathogenic
print("Bayesian Points:", report.bayesian_points)      # 14
print("Triggered Rules:", report.categorical_triggered_rules)

# 2. Splice Impact Prediction
splice_pred = SpliceImpactPredictor().predict("VAR1", position=1, ref_base="G", alt_base="A", splice_distance=1)
print("Splice Impact:", splice_pred.predicted_effect, splice_pred.evidence_level)

# 3. Pharmacogenomics Annotation
pgx_annot = PharmacoGenomicsAnnotator().annotate("CYP2C19", "c.681G>A")
print(f"CPIC Level: {pgx_annot.cpic_level} for {pgx_annot.drug_names}")
```

---

## Running Tests

Run the test suite using standard `unittest` or `pytest`:

```bash
python -m unittest discover tests
# or
pytest -v
```

