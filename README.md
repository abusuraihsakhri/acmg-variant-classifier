# Acmg Variant Classifier

> **Domain:** Clinical Decision Support & Biomedical Computing  
> **Reference Guidelines & Standards:** `Standard Clinical Formulations & ISO/IEC Quality Frameworks`

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## 📖 What It Does

ACMG/AMP 2015 & ClinGen / Tavtigian 2020 Variant Classifier & Functional Genomics Annotation.

Provides:
- Comprehensive 28-code ACMG/AMP 2015 pathogenicity evaluation.
- Tavtigian et al. (2018/2020) Bayesian point scoring & posterior probability model.
- Automated population allele frequency checking (gnomAD / 1000 Genomes thresholds).
- Splice site impact prediction & functional domain hotspot mapping.
- CPIC pharmacogenomic variant annotation.

ACMG/AMP Variant Classifier: Splice Variant Impact Predictor & Pharmaco-Genomics Annotation.
Implements in-silico splice prediction, PGx CPIC level annotation, and phenotype-variant correlation.

---

## ⚙️ Key Capabilities & Algorithmic Modules

### 🔬 Core Algorithmic & Evaluation Engines

- **`SplicePrediction`** — dedicated module for splice prediction evaluation and state verification.
- **`PGxAnnotation`** — dedicated module for p gx annotation evaluation and state verification.
- **`PhenotypeCorrelation`** — dedicated module for phenotype correlation evaluation and state verification.
- **`FunctionalDomainHit`** — dedicated module for functional domain hit evaluation and state verification.
- **`SpliceImpactPredictor`**: In-silico splice variant impact prediction using position-weight matrices.
- **`PharmacoGenomicsAnnotator`**: ClinGen/CPIC pharmacogenomics annotation for variants.

---

## 📐 Mathematical Formulation & Logic

```text
  delta_score = 0.0
  delta_score = 0.9 if ref_base in ("G", "A", "T") else 0.7
  delta_score = 0.5
  delta_score = 0.2
  delta_score = 0.3
```

---

## 💻 CLI Quickstart & Usage

### 1. Guided Interactive Mode
```bash
python cli.py
```

### 2. Direct Parameterized Evaluation
```bash
python cli.py --evidence <value> --input <value> --vcf <value> --interactive <value>
```

### Parameter Reference
- `--evidence`: Specifies input measurement or parameter value.
- `--input`: Specifies input measurement or parameter value.
- `--vcf`: Specifies input measurement or parameter value.
- `--interactive`: Specifies input measurement or parameter value.
- `--splice`: Specifies input measurement or parameter value.
- `--pgx`: Specifies input measurement or parameter value.
- `--domain`: Specifies input measurement or parameter value.
- `--variant-id`: Specifies input measurement or parameter value.
- `--gene`: Specifies input measurement or parameter value.
- `--af`: Specifies input measurement or parameter value.

### Input Data Schema

| Field | Description | Requirement |
|:------|:------------|:------------|
| `suite_name` | Parameter / observation metric | Required |
| `system_slug` | Parameter / observation metric | Required |
| `standard_reference` | Parameter / observation metric | Required |
| `benchmark_variants` | Parameter / observation metric | Required |

---

## 🛡️ Security & Enterprise Architecture

* **Zero-PHI Outbound Interceptor:** Active AST and regex inspection blocking SSNs, MRNs, phone numbers, and patient identifiers.
* **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation and state transition.
* **Air-Gapped LLM Reasoning Adapter:** Agnostic integration for local Ollama instances (`llama3`, `mistral`), Claude 3.5 Sonnet, GPT-4o, and deterministic test mocks.
* **Active Learning Bayesian Calibration:** Dynamic tracker updating worker reliability weights and monitoring Brier calibration drift.
* **FastAPI & Prometheus Telemetry:** Exposes OpenAPI 3.1 REST endpoints and operational Prometheus metrics (`/metrics`).

---

## 🧪 Testing & Verification

Run the automated test suite:

```bash
pytest -v
```

Execute high-throughput batch simulation benchmarks:

```bash
python simulator.py --tasks 1000 --concurrency 8
```

---

## 🐳 Container Deployment

```bash
docker build -t acmg-variant-classifier .
docker run -p 8000:8000 acmg-variant-classifier
```
