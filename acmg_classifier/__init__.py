"""ACMG Variant Classifier.

Implements the 2015 ACMG/AMP sequence variant interpretation guidelines
(Richards et al., Genet Med 2015) combining-rule table, together with the
Tavtigian et al. 2020 (Genet Med) Bayesian point-scoring refinement, to
classify germline variants into the five-tier system: Pathogenic, Likely
Pathogenic, Uncertain Significance, Likely Benign, Benign.
"""

__version__ = "1.0.0"

from .classifier import classify_variant, ClassificationReport
from .bayesian import bayesian_classification, posterior_probability, BAYESIAN_POINTS
from .evidence import validate_codes, count_by_category, VALID_CODES, EVIDENCE_CODES
from .frequency import check_frequency, DEFAULT_BA1_THRESHOLD, DEFAULT_PM2_THRESHOLD
from .rules import categorical_classification
from .report import report_to_dict, reports_to_json, report_to_text, reports_to_text, reports_to_csv
from .security import safe_resolve_path, sanitize_csv_cell, validate_numeric_range, WINDOWS_RESERVED_NAMES

__all__ = [
    "__version__",
    "classify_variant",
    "ClassificationReport",
    "bayesian_classification",
    "posterior_probability",
    "BAYESIAN_POINTS",
    "validate_codes",
    "count_by_category",
    "VALID_CODES",
    "EVIDENCE_CODES",
    "check_frequency",
    "DEFAULT_BA1_THRESHOLD",
    "DEFAULT_PM2_THRESHOLD",
    "categorical_classification",
    "report_to_dict",
    "reports_to_json",
    "report_to_text",
    "reports_to_text",
    "reports_to_csv",
    "safe_resolve_path",
    "sanitize_csv_cell",
    "validate_numeric_range",
    "WINDOWS_RESERVED_NAMES",
]
