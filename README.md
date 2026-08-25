# ACMG/AMP 2015 & Tavtigian 2020 Variant Classifier

Production-grade clinical genomics variant pathogenicity classifier and functional annotation engine implementing the **ACMG/AMP 2015 standards** (*Richards et al., Genetics in Medicine*) and the **ClinGen / Tavtigian 2020 Bayesian point-scoring refinement** (*Tavtigian et al., Human Mutation*).

---

## Key Capabilities

1. **Full 28 ACMG/AMP Evidence Criteria**:
   - **Pathogenic Very Strong**: `PVS1` (null variants, frameshifts, nonsense, canonical splice disruption).
   - **Pathogenic Strong**: `PS1`, `PS2`, `PS3`, `PS4` (well-established functional assays, de novo, identical amino acid changes).
   - **Pathogenic Moderate**: `PM1`, `PM2`, `PM3`, `PM4`, `PM5`, `PM6` (mutational hotspots, population absence, in-trans segregation).
   - **Pathogenic Supporting**: `PP1`, `PP2`, `PP3`, `PP4`, `PP5` (in-silico predictions, phenotype specificity).
   - **Benign Stand-Alone**: `BA1` (allele frequency $\ge 5\%$).
   - **Benign Strong**: `BS1`, `BS2`, `BS3`, `BS4` (disease allele frequency threshold, healthy adult controls, functional assays).
   - **Benign Supporting**: `BP1`, `BP2`, `BP3`, `BP4`, `BP5`, `BP6`, `BP7` (in-silico neutral, synonymous with no splice effect).

2. **Dual-Classification Engine**:
   - **Categorical Rule Table**: Standard Richards et al. 2015 multi-tier logic matrices.
   - **Bayesian Odds & Point Scaling**: Tavtigian et al. exponential scoring (+8 for Very Strong, +4 for Strong, +2 for Moderate, +1 for Supporting; -1 for Benign Supporting, -4 for Benign Strong).

3. **Population Allele Frequency Cross-Checking**:
   - Automatic evaluation against gnomAD and 1000 Genomes frequency thresholds.
   - Detection of contradiction between manual inputs (e.g. asserting `PM2` on a common polymorphism).

4. **In-Silico Splice & Functional Genomic Modules**:
   - Splice junction disruption score and delta calculation.
   - CPIC (Clinical Pharmacogenetics Implementation Consortium) level A/B annotation.
   - Protein functional domain and hotspot mapping (BRCA1, TP53, CFTR, BRAF).

---

## Mathematical & Bayesian Formulation

The Tavtigian Bayesian framework models evidence combinations as multiplicative odds:

$$\text{Odds}_{\text{Path}} = \frac{P(\text{Pathogenic} \mid E)}{1 - P(\text{Pathogenic} \mid E)} = \text{Odds}_{\text{prior}} \times \prod_{i} (\text{OddsPath}_i)^{w_i}$$

Where individual evidence strengths map to integer points:
- **Very Strong ($w=+8$)**: $\text{OddsPath} \approx 350:1$
- **Strong ($w=+4$)**: $\text{OddsPath} \approx 18.7:1$
- **Moderate ($w=+2$)**: $\text{OddsPath} \approx 4.33:1$
- **Supporting ($w=+1$)**: $\text{OddsPath} \approx 2.08:1$
- **Benign Supporting ($w=-1$)**: $\text{OddsPath} \approx 0.48:1$
- **Benign Strong ($w=-4$)**: $\text{OddsPath} \approx 0.053:1$

Points scale into clinical tiers:
- $\ge 10\text{ points} \implies \text{Pathogenic}$ ($P \ge 0.99$)
- $6 \text{ to } 9\text{ points} \implies \text{Likely Pathogenic}$ ($0.90 \le P < 0.99$)
- $0 \text{ to } 5\text{ points} \implies \text{Uncertain Significance (VUS)}$
- $-1 \text{ to } -6\text{ points} \implies \text{Likely Benign}$ ($0.01 < P \le 0.10$)
- $\le -7\text{ points} \text{ or } \text{BA1} \implies \text{Benign}$ ($P \le 0.01$)

---

## Command Line Interface (CLI)

### 1. Single Variant Classification
```bash
python cli.py -e PVS1 PS3 PM2 --variant-id "BRCA1:c.5266dupC" --af 0.00002
```

### 2. Output as JSON
```bash
python cli.py -e PVS1 PM1 PP3 --variant-id "TP53:c.524G>A" --format json
```

### 3. Batch Processing (TSV / CSV)
```bash
python cli.py -i variants.tsv -o results.json --format json
```

### 4. VCF File Annotation
```bash
python cli.py --vcf sample.vcf -o vcf_classified.txt
```

### 5. In-Silico Splice Prediction
```bash
python cli.py --splice --variant-id "MSH2:c.942+1G>A" --pos 1 --ref G --alt A --region-type intron --distance 1
```

### 6. Pharmacogenomics (CPIC) Annotation
```bash
python cli.py --pgx --gene CYP2C19 --variant-id "c.681G>A" --variant-type nonsense
```

### 7. Interactive Mode
```bash
python cli.py --interactive
```

---

## Running the Unit Test Suite

Execute the unit test suite with 100% standard library Python:

```bash
python -m unittest discover -s tests -p "test_*.py" -v
# or
python test_classifier.py
```
