<div align="center">

# 📦 ACMG Variant Classifier & Functional Genomics Annotator

**ACMG/AMP germline evidence codes and population allele frequencies evaluated via Richards et al. combination rules and Tavtigian Bayesian point modeling to produce five-tier pathogenicity classifications, splice disruption scores, and CPIC pharmacogenomic annotations in text, JSON, or CSV with zero external runtime dependencies.**

[![Python Version](https://img.shields.io/badge/python-3.10%2B-blue.svg)]()
[![Zero Dependencies](https://img.shields.io/badge/dependencies-0%20(stdlib%20only)-success.svg)]()
[![Test Status](https://img.shields.io/badge/tests-100%25%20passed-brightgreen.svg)]()
[![License](https://img.shields.io/badge/license-MIT-purple.svg)](LICENSE)

</div>

---

## 📖 Executive Summary

Clinical genetics pipelines require consistent, reproducible interpretation of DNA sequence variants across diagnostic and research workflows. This library evaluates germline variants using the 28 criteria of the 2015 ACMG/AMP framework alongside the 2020 Tavtigian Bayesian point-scoring refinement, resolving conflicting evidence tiers into calibrated posterior probabilities. All algorithms, population frequency cross-checks, splice disruption models, and CPIC pharmacogenomic mappings run entirely within the Python standard library with zero third-party runtime dependencies.

---

## 🧬 Architecture & Logic Flowchart

```
[Variant Input (CLI / TSV / CSV / VCF)]
               │
               ▼
[Defensive Input Validation]
  ├── Path Resolution & Directory Confinement (CWE-22)
  ├── Windows Reserved Device Neutralization
  └── Strict Numeric Range & Boundary Checks (NaN / Inf / bool)
               │
               ▼
[Allele Frequency Cross-Check]
  ├── BA1 Stand-Alone Benign Check (AF > 0.05)
  ├── PM2 Absent/Rare Control Check (AF <= 0.0001)
  └── Contradiction & Discordance Detection
               │
               ▼
[Dual Classification Engines]
  ├── ACMG/AMP 2015 Combining Rules (Table 5)
  └── Tavtigian 2020 Bayesian Point Sum & Posterior Probability
               │
               ▼
[Genomic & Clinical Annotation]
  ├── Splice Impact Predictor (Position-Weight Donor/Acceptor)
  ├── Functional Domain & Mutational Hotspot Mapper
  └── CPIC Pharmacogenomics Annotator (Level A/B Actionable Drugs)
               │
               ▼
[Sanitized Multi-Format Output]
  ├── Human-Readable Formatted Terminal Text
  ├── Strict RFC 8259 JSON (NaN / Inf neutralized)
  └── Formula-Sanitized CSV (CWE-1236 Neutralization)
```

---

## 🚀 Installation

```bash
git clone https://github.com/abusuraihsakhri/acmg-variant-classifier.git
cd acmg-variant-classifier
pip install -e .
```

For development and running test suites:
```bash
pip install -r requirements-dev.txt
```

---

## 💻 CLI Usage & Command Reference

### 1. Single Variant Evaluation
Evaluate evidence criteria for a single variant:
```bash
python cli.py --evidence PVS1 PS3 PM2 --variant-id "BRCA1:c.5266dupC"
```
*Output:*
```text
Variant: BRCA1:c.5266dupC
Final classification: Pathogenic
  determined by: ACMG/AMP 2015 categorical rule table
Categorical (ACMG/AMP 2015): Pathogenic
Bayesian (Tavtigian 2020): Pathogenic (14 points)
Input evidence codes: PVS1, PS3, PM2
Evidence counts: PVS=1, PS=1, PM=1
Triggered categorical rules:
  - Pathogenic (i-a): 1 PVS1 AND >=1 PS
  - Likely Pathogenic (1): 1 PVS1 AND 1 PM
  - Likely Pathogenic (2): 1 PS AND 1-2 PM
Rationale: Evidence tallied: PVS=1, PS=1, PM=1. Categorical rules (Richards et al. 2015) triggered: Pathogenic (i-a): 1 PVS1 AND >=1 PS; Likely Pathogenic (1): 1 PVS1 AND 1 PM; Likely Pathogenic (2): 1 PS AND 1-2 PM. Tavtigian 2020 Bayesian point total: 14 points -> Pathogenic. Final classification: Pathogenic (determined by ACMG/AMP 2015 categorical rule table).
```

*Output as JSON:*
```bash
python cli.py --evidence PVS1 PS3 PM2 --variant-id "BRCA1:c.5266dupC" --json
```
```json
[
  {
    "variant_id": "BRCA1:c.5266dupC",
    "input_codes": ["PVS1", "PS3", "PM2"],
    "unknown_codes": [],
    "effective_codes": ["PM2", "PS3", "PVS1"],
    "auto_added_codes": [],
    "evidence_counts": {"PVS": 1, "PS": 1, "PM": 1, "PP": 0, "BA": 0, "BS": 0, "BP": 0},
    "categorical_label": "Pathogenic",
    "categorical_triggered_rules": [
      "Pathogenic (i-a): 1 PVS1 AND >=1 PS",
      "Likely Pathogenic (1): 1 PVS1 AND 1 PM",
      "Likely Pathogenic (2): 1 PS AND 1-2 PM"
    ],
    "categorical_conflict": false,
    "bayesian_points": 14,
    "bayesian_label": "Pathogenic",
    "final_classification": "Pathogenic",
    "final_source": "ACMG/AMP 2015 categorical rule table",
    "warnings": [],
    "rationale": "Evidence tallied: PVS=1, PS=1, PM=1..."
  }
]
```

### 2. Batch Tabular Evaluation (CSV / TSV)
Process tabular variant files with automatic delimiter detection and per-row crash isolation:
```bash
python cli.py -i sample_variants.csv --format csv
```
*Output:*
```csv
variant_id,final_classification,final_source,categorical_label,bayesian_points,bayesian_label,input_codes,effective_codes,auto_added_codes,warnings,rationale
VAR_BRCA1_01,Pathogenic,ACMG/AMP 2015 categorical rule table,Pathogenic,14,Pathogenic,PVS1; PS3; PM2,PM2; PS3; PVS1,,,"Evidence tallied: PVS=1, PS=1, PM=1..."
VAR_TP53_01,Pathogenic,ACMG/AMP 2015 categorical rule table,Pathogenic,13,Pathogenic,PS1; PS3; PM1; PM2; PP3,PM1; PM2; PP3; PS1; PS3,,,"Evidence tallied: PS=2, PM=2, PP=1..."
VAR_COMMON_01,Benign,ACMG/AMP 2015 categorical rule table,Benign,0,Benign,BA1,BA1,,,Evidence tallied: BA=1...
```

### 3. Splice Site Disruption Prediction
Predict splice donor/acceptor alterations using position-weight matrices:
```bash
python cli.py --splice --pos 1 --ref G --alt A --distance 1 --json
```
```json
{
  "variant_id": "VAR-001",
  "splice_site_type": "canonical_acceptor",
  "predicted_effect": "Loss of canonical splice site",
  "delta_score": 0.9,
  "conservation_score": 0.95,
  "evidence_level": "strong"
}
```

### 4. Pharmacogenomics (CPIC) Annotation
Annotate actionable clinical guidelines across CPIC Level A/B genes:
```bash
python cli.py --pgx --gene CYP2C19 --variant-id "c.681G>A" --json
```
```json
{
  "variant_id": "c.681G>A",
  "gene": "CYP2C19",
  "cpic_level": "A",
  "drug_names": ["clopidogrel", "omeprazole", "voriconazole"],
  "phenotype": "CYP2C19 poor metabolizer",
  "recommendation": "May affect CYP2C19 function. Consider phenotyping if clinically relevant.",
  "evidence_url": "https://www.pharmgkb.org/gene/CYP2C19/clinicalAnnotation"
}
```

### 5. Functional Domain & Hotspot Mapping
Cross-reference variant amino acid coordinates against domain boundaries:
```bash
python cli.py --domain --gene TP53 --pos 175 --json
```
```json
[
  {
    "domain_name": "DNA-binding domain",
    "residue_range": [102, 292],
    "variant_position": 175,
    "in_domain": true,
    "hot_spot": true
  }
]
```

### 6. Interactive Terminal Wizard
```bash
python cli.py --interactive
```

---

## 🐍 Python API Reference

```python
from acmg_classifier import classify_variant, reports_to_json, reports_to_csv
from variant_analysis import SpliceImpactPredictor, PharmacoGenomicsAnnotator, FunctionalDomainMapper

# 1. Classify a clinical variant
report = classify_variant(
    codes=["PVS1", "PS3", "PM2"],
    variant_id="BRCA1:c.5266dupC",
    population_af=0.00002
)
print("Classification:", report.final_classification)  # Pathogenic
print("Bayesian Points:", report.bayesian_points)      # 14
print("Triggered Rules:", report.categorical_triggered_rules)

# 2. Predict splice junction disruption
splice_result = SpliceImpactPredictor().predict(
    variant_id="SP_01",
    position=1,
    ref_base="G",
    alt_base="A",
    splice_distance=1
)
print("Splice Impact:", splice_result.predicted_effect, splice_result.evidence_level)

# 3. Retrieve CPIC pharmacogenomics guidance
pgx_result = PharmacoGenomicsAnnotator().annotate("CYP2D6", "CYP2D6*4", variant_type="nonsense")
print(f"CPIC Level: {pgx_result.cpic_level} for {pgx_result.drug_names}")
```

---

## 📊 Reference Tables & Specifications

### ACMG/AMP 2015 & Tavtigian 2020 Point Equivalencies

| Strength Category | Code Prefix | Tavtigian Points | Direction | Description |
|:------------------|:------------|:-----------------|:----------|:------------|
| Very Strong       | `PVS1`      | +8               | Pathogenic| Null variant where loss of function is known disease mechanism |
| Strong            | `PS1`-`PS4` | +4               | Pathogenic| Established functional studies, de novo, amino acid match |
| Moderate          | `PM1`-`PM6` | +2               | Pathogenic| Mutational hot spot, absent from controls, protein length change |
| Supporting        | `PP1`-`PP5` | +1               | Pathogenic| Computational in-silico evidence, co-segregation, phenotype match |
| Stand-Alone Benign| `BA1`       | N/A (Decisive)   | Benign    | Population allele frequency > 5% (`> 0.05`) |
| Strong Benign     | `BS1`-`BS4` | -4               | Benign    | Greater frequency than expected, adult healthy controls |
| Supporting Benign | `BP1`-`BP7` | -1               | Benign    | Computational benign, synonymous non-conserved, alternative cause |

### Bayesian Classification Thresholds

| Classification Tier | Bayesian Point Range | Posterior Probability Range |
|:--------------------|:---------------------|:----------------------------|
| Pathogenic          | $\ge +10$            | $\ge 0.99$                  |
| Likely Pathogenic   | $+6$ to $+9$         | $0.90$ to $0.98$            |
| Uncertain (VUS)     | $0$ to $+5$          | $0.10$ to $0.89$            |
| Likely Benign       | $-1$ to $-6$         | $0.01$ to $0.09$            |
| Benign              | $\le -7$ or `BA1`    | $< 0.01$                    |

---

## 🛡️ Security & Defensive Architecture

| Threat Category | Risk | Implemented Mitigation |
|:----------------|:-----|:-----------------------|
| CSV Injection (CWE-1236) | Malicious formulas (`=`, `+`, `@`, `\t`, `\r`, `-cmd`) executing in Excel | `sanitize_csv_cell` prepends single quote (`'`) to dangerous characters while preserving numbers |
| Type & Boundary Bypass (CWE-1287) | `NaN`, `Inf`, and `bool` bypassing inequality thresholds | Strict `validate_numeric_range` enforcing non-boolean, finite numbers within bounds |
| Directory Traversal & Device Hang (CWE-22) | Escape outside workdir (`../`) and Windows device lockups (`CON`, `PRN`, `NUL`) | `safe_resolve_path` with `os.path.commonpath` confinement and device name blacklisting |
| Stack Trace Leakage (CWE-209) | Raw tracebacks dumping server paths and internal structures | Global CLI exception handling returning sanitized exit codes and machine-readable JSON errors |
| Batch Cascade Failure (A08:2025) | Single corrupted row aborting large-scale variant batches | Per-row try/except isolation recording row-level error reports without pipeline halting |
| Non-Compliant JSON (RFC 8259) | Unquoted `NaN` or `Infinity` tokens breaking strict parsers | `allow_nan=False` serialization and finite float sanitization across all endpoints |

---

## 🧪 Verification & Testing

Execute the test suite covering functional logic, edge cases, and defensive boundaries:

```bash
pytest -v
```

Execute standalone test runner:
```bash
python test_classifier.py
```

Run CLI smoke verification:
```bash
python cli.py --evidence PVS1 PS3 PM2 --json
python cli.py -i sample_variants.csv --format csv
```

---

## 📄 License

MIT License. See [LICENSE](LICENSE) for details.
