"""
ACMG/AMP 2015 & ClinGen / Tavtigian 2020 Variant Classifier & Functional Genomics Annotation.

Provides:
- Comprehensive 28-code ACMG/AMP 2015 pathogenicity evaluation.
- Tavtigian et al. (2018/2020) Bayesian point scoring & posterior probability model.
- Automated population allele frequency checking (gnomAD / 1000 Genomes thresholds).
- Splice site impact prediction & functional domain hotspot mapping.
- CPIC pharmacogenomic variant annotation.
"""

from acmg_classifier.evidence import (
    VALID_CODES,
    CATEGORIES_PATHOGENIC,
    CATEGORIES_BENIGN,
    EVIDENCE_CODES,
    validate_codes,
    count_by_category,
)
from acmg_classifier.rules import categorical_classification
from acmg_classifier.bayesian import (
    bayesian_classification,
    posterior_probability,
    BAYESIAN_POINTS,
)
from acmg_classifier.frequency import (
    check_frequency,
    DEFAULT_BA1_THRESHOLD,
    DEFAULT_PM2_THRESHOLD,
)
from acmg_classifier.classifier import (
    classify_variant,
    ClassificationReport,
)
from acmg_classifier.report import (
    report_to_dict,
    reports_to_json,
    report_to_text,
    reports_to_text,
)
from variant_analysis import (
    SplicePrediction,
    SpliceImpactPredictor,
    PGxAnnotation,
    PharmacoGenomicsAnnotator,
    FunctionalDomainHit,
    FunctionalDomainMapper,
    PhenotypeCorrelation,
)

__all__ = [
    "VALID_CODES",
    "CATEGORIES_PATHOGENIC",
    "CATEGORIES_BENIGN",
    "EVIDENCE_CODES",
    "validate_codes",
    "count_by_category",
    "categorical_classification",
    "bayesian_classification",
    "posterior_probability",
    "BAYESIAN_POINTS",
    "check_frequency",
    "DEFAULT_BA1_THRESHOLD",
    "DEFAULT_PM2_THRESHOLD",
    "classify_variant",
    "ClassificationReport",
    "report_to_dict",
    "reports_to_json",
    "report_to_text",
    "reports_to_text",
    "SplicePrediction",
    "SpliceImpactPredictor",
    "PGxAnnotation",
    "PharmacoGenomicsAnnotator",
    "FunctionalDomainHit",
    "FunctionalDomainMapper",
    "PhenotypeCorrelation",
]
