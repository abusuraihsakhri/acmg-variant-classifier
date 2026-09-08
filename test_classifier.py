"""
ACMG Variant Classifier & Annotation Suite Tests.
Runs all unit test cases for categorical rules, Bayesian points, allele frequency filtering,
splice disruption prediction, and CPIC pharmacogenomics.
"""

import sys
import unittest
from pathlib import Path

# Ensure package path is resolved
sys.path.insert(0, str(Path(__file__).parent))

from tests.test_acmg_variant_classifier import (
    TestEvidenceValidation,
    TestCategoricalPathogenicRules,
    TestCategoricalLikelyPathogenicRules,
    TestCategoricalBenignAndVUS,
    TestBayesianScoring,
    TestAlleleFrequencyRules,
    TestSpliceAndPGxAndDomain,
    TestCLIAndBatch,
    TestInputValidationAndSecurity,
    TestDefensiveHardeningAndSecurity,
)

if __name__ == "__main__":
    unittest.main()
