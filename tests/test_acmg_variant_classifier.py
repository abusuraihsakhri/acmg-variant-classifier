"""
Comprehensive Unit Test Suite for ACMG/AMP 2015 & Tavtigian 2020 Variant Classifier.
Covers categorical rules, Bayesian point modeling, allele frequency constraints,
splice impact predictions, CPIC pharmacogenomics, functional domain mapping, and CLI interfaces.
"""

import io
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

# Add project root to sys.path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from acmg_classifier import (
    evidence,
    rules,
    bayesian,
    frequency,
    classifier,
    report as report_mod,
)
from acmg_classifier.classifier import classify_variant
from acmg_variant_classifier import (
    SpliceImpactPredictor,
    PharmacoGenomicsAnnotator,
    FunctionalDomainMapper,
)
import cli


class TestEvidenceValidation(unittest.TestCase):
    def test_valid_codes_normalized(self):
        codes, unknown = evidence.validate_codes(["ps1", " PM2 ", "PVS1", "bp4"])
        self.assertEqual(codes, ["PS1", "PM2", "PVS1", "BP4"])
        self.assertEqual(unknown, [])

    def test_unknown_code_detected(self):
        codes, unknown = evidence.validate_codes(["PS1", "PZ9", "INVALID"])
        self.assertEqual(unknown, ["PZ9", "INVALID"])

    def test_duplicate_codes_deduped(self):
        codes, unknown = evidence.validate_codes(["PS1", "PS1", "PM2", "PM2"])
        self.assertEqual(codes, ["PS1", "PM2"])

    def test_unknown_code_raises_in_classifier(self):
        with self.assertRaises(ValueError):
            classify_variant(["PS1", "NOT_A_CODE"])

    def test_empty_codes_handling(self):
        codes, unknown = evidence.validate_codes([])
        self.assertEqual(codes, [])
        self.assertEqual(unknown, [])


class TestCategoricalPathogenicRules(unittest.TestCase):
    def test_pathogenic_pvs1_plus_1ps(self):
        counts = evidence.count_by_category(["PVS1", "PS3"])
        res = rules.categorical_classification(counts)
        self.assertEqual(res["label"], "Pathogenic")
        self.assertTrue(res["pathogenic_tier_hit"])

    def test_pathogenic_pvs1_plus_2pm(self):
        counts = evidence.count_by_category(["PVS1", "PM1", "PM2"])
        res = rules.categorical_classification(counts)
        self.assertEqual(res["label"], "Pathogenic")

    def test_pathogenic_pvs1_plus_1pm_1pp(self):
        counts = evidence.count_by_category(["PVS1", "PM1", "PP1"])
        res = rules.categorical_classification(counts)
        self.assertEqual(res["label"], "Pathogenic")

    def test_pathogenic_pvs1_plus_2pp(self):
        counts = evidence.count_by_category(["PVS1", "PP1", "PP2"])
        res = rules.categorical_classification(counts)
        self.assertEqual(res["label"], "Pathogenic")

    def test_pathogenic_2ps(self):
        counts = evidence.count_by_category(["PS1", "PS2"])
        res = rules.categorical_classification(counts)
        self.assertEqual(res["label"], "Pathogenic")

    def test_pathogenic_1ps_plus_3pm(self):
        counts = evidence.count_by_category(["PS1", "PM1", "PM2", "PM3"])
        res = rules.categorical_classification(counts)
        self.assertEqual(res["label"], "Pathogenic")

    def test_pathogenic_1ps_plus_2pm_plus_2pp(self):
        counts = evidence.count_by_category(["PS1", "PM1", "PM2", "PP1", "PP2"])
        res = rules.categorical_classification(counts)
        self.assertEqual(res["label"], "Pathogenic")

    def test_pathogenic_1ps_plus_1pm_plus_4pp(self):
        counts = evidence.count_by_category(["PS1", "PM1", "PP1", "PP2", "PP3", "PP4"])
        res = rules.categorical_classification(counts)
        self.assertEqual(res["label"], "Pathogenic")


class TestCategoricalLikelyPathogenicRules(unittest.TestCase):
    def test_likely_pathogenic_1pvs1_1pm(self):
        counts = evidence.count_by_category(["PVS1", "PM1"])
        res = rules.categorical_classification(counts)
        self.assertEqual(res["label"], "Likely Pathogenic")

    def test_likely_pathogenic_1ps_1pm(self):
        counts = evidence.count_by_category(["PS1", "PM1"])
        res = rules.categorical_classification(counts)
        self.assertEqual(res["label"], "Likely Pathogenic")

    def test_likely_pathogenic_1ps_2pp(self):
        counts = evidence.count_by_category(["PS1", "PP1", "PP2"])
        res = rules.categorical_classification(counts)
        self.assertEqual(res["label"], "Likely Pathogenic")

    def test_likely_pathogenic_3pm(self):
        counts = evidence.count_by_category(["PM1", "PM2", "PM3"])
        res = rules.categorical_classification(counts)
        self.assertEqual(res["label"], "Likely Pathogenic")

    def test_likely_pathogenic_2pm_2pp(self):
        counts = evidence.count_by_category(["PM1", "PM2", "PP1", "PP2"])
        res = rules.categorical_classification(counts)
        self.assertEqual(res["label"], "Likely Pathogenic")

    def test_likely_pathogenic_1pm_4pp(self):
        counts = evidence.count_by_category(["PM1", "PP1", "PP2", "PP3", "PP4"])
        res = rules.categorical_classification(counts)
        self.assertEqual(res["label"], "Likely Pathogenic")


class TestCategoricalBenignAndVUS(unittest.TestCase):
    def test_benign_ba1_standalone(self):
        counts = evidence.count_by_category(["BA1"])
        res = rules.categorical_classification(counts)
        self.assertEqual(res["label"], "Benign")

    def test_benign_2bs(self):
        counts = evidence.count_by_category(["BS1", "BS2"])
        res = rules.categorical_classification(counts)
        self.assertEqual(res["label"], "Benign")

    def test_likely_benign_1bs_1bp(self):
        counts = evidence.count_by_category(["BS1", "BP1"])
        res = rules.categorical_classification(counts)
        self.assertEqual(res["label"], "Likely Benign")

    def test_likely_benign_2bp(self):
        counts = evidence.count_by_category(["BP1", "BP2"])
        res = rules.categorical_classification(counts)
        self.assertEqual(res["label"], "Likely Benign")

    def test_vus_insufficient_evidence(self):
        counts = evidence.count_by_category(["PP3"])
        res = rules.categorical_classification(counts)
        self.assertEqual(res["label"], "Uncertain Significance")

    def test_conflict_pathogenic_and_benign(self):
        counts = evidence.count_by_category(["PVS1", "PS1", "BS1", "BS2"])
        res = rules.categorical_classification(counts)
        self.assertTrue(res["conflict"])


class TestBayesianScoring(unittest.TestCase):
    def test_bayesian_pathogenic_threshold(self):
        counts = evidence.count_by_category(["PVS1", "PS1"])  # 8 + 4 = 12 >= 10
        res = bayesian.bayesian_classification(counts, ba1_present=False)
        self.assertEqual(res["points"], 12)
        self.assertEqual(res["label"], "Pathogenic")

    def test_bayesian_likely_pathogenic_threshold(self):
        counts = evidence.count_by_category(["PS1", "PM1"])  # 4 + 2 = 6 (6-9 -> LP)
        res = bayesian.bayesian_classification(counts, ba1_present=False)
        self.assertEqual(res["points"], 6)
        self.assertEqual(res["label"], "Likely Pathogenic")

    def test_bayesian_benign_threshold(self):
        counts = evidence.count_by_category(["BS1", "BS2"])  # -4 + -4 = -8 (<= -7 -> Benign)
        res = bayesian.bayesian_classification(counts, ba1_present=False)
        self.assertLessEqual(res["points"], -7)
        self.assertEqual(res["label"], "Benign")

    def test_bayesian_ba1_stand_alone(self):
        counts = evidence.count_by_category(["BA1"])
        res = bayesian.bayesian_classification(counts, ba1_present=True)
        self.assertEqual(res["label"], "Benign")

    def test_posterior_probability_calculation(self):
        prob_path = bayesian.posterior_probability(points=10, prior=0.10)
        self.assertGreater(prob_path, 0.95)

        prob_benign = bayesian.posterior_probability(points=-7, prior=0.10)
        self.assertLess(prob_benign, 0.05)


class TestAlleleFrequencyRules(unittest.TestCase):
    def test_auto_ba1_when_af_exceeds_threshold(self):
        effective_codes, warnings, auto_added = frequency.check_frequency(
            ["PP1"], population_af=0.06, ba1_threshold=0.05, pm2_threshold=0.0001, auto_apply=True
        )
        self.assertIn("BA1", effective_codes)
        self.assertIn("BA1", auto_added)

    def test_auto_pm2_when_af_below_threshold(self):
        effective_codes, warnings, auto_added = frequency.check_frequency(
            ["PP1"], population_af=0.00002, ba1_threshold=0.05, pm2_threshold=0.0001, auto_apply=True
        )
        self.assertIn("PM2", effective_codes)
        self.assertIn("PM2", auto_added)

    def test_frequency_contradiction_warning(self):
        effective_codes, warnings, auto_added = frequency.check_frequency(
            ["PM2"], population_af=0.08, ba1_threshold=0.05, pm2_threshold=0.0001, auto_apply=False
        )
        self.assertTrue(any("contradicts" in w or "BA1" in w for w in warnings))


class TestSpliceAndPGxAndDomain(unittest.TestCase):
    def test_splice_canonical_donor_impact(self):
        predictor = SpliceImpactPredictor()
        res = predictor.predict(
            variant_id="V_SPLICE_1",
            position=1,
            ref_base="G",
            alt_base="A",
            intron_exon="intron",
            splice_distance=1,
        )
        self.assertEqual(res.splice_site_type, "canonical_acceptor")
        self.assertGreaterEqual(res.delta_score, 0.7)
        self.assertEqual(res.evidence_level, "strong")

    def test_splice_deep_intronic(self):
        predictor = SpliceImpactPredictor()
        res = predictor.predict(
            variant_id="V_DEEP_INTRON",
            position=100,
            ref_base="G",
            alt_base="C",
            intron_exon="intron",
            splice_distance=25,
        )
        self.assertEqual(res.splice_site_type, "deep_intronic")
        self.assertEqual(res.evidence_level, "supporting")

    def test_pgx_cyp2d6_annotation(self):
        annotator = PharmacoGenomicsAnnotator()
        res = annotator.annotate("CYP2D6", "CYP2D6*4", "nonsense")
        self.assertEqual(res.cpic_level, "A")
        self.assertIn("codeine", res.drug_names)
        self.assertIn("CYP2D6", res.recommendation)

    def test_pgx_unknown_gene(self):
        annotator = PharmacoGenomicsAnnotator()
        res = annotator.annotate("UNKNOWN_GENE", "VAR_X", "missense")
        self.assertEqual(res.cpic_level, "No Data")

    def test_domain_mapping_brca1_ring(self):
        mapper = FunctionalDomainMapper()
        hits = mapper.map_variant("BRCA1", "BRCA1:p.C64G", position=64)
        self.assertTrue(any(h.domain_name == "RING domain" and h.in_domain for h in hits))

    def test_domain_mapping_tp53_hotspot(self):
        mapper = FunctionalDomainMapper()
        hits = mapper.map_variant("TP53", "TP53:p.R175H", position=175)
        self.assertTrue(any(h.hot_spot for h in hits))


class TestCLIAndBatch(unittest.TestCase):
    def test_cli_single_variant_json(self):
        out = io.StringIO()
        err = io.StringIO()
        old_stdout = sys.stdout
        sys.stdout = out
        try:
            exit_code = cli.main(["-e", "PVS1", "PS3", "--format", "json", "--variant-id", "TEST-01"])
            self.assertEqual(exit_code, 0)
        finally:
            sys.stdout = old_stdout

        payload = json.loads(out.getvalue())
        self.assertEqual(payload[0]["final_classification"], "Pathogenic")
        self.assertEqual(payload[0]["variant_id"], "TEST-01")

    def test_cli_tsv_batch(self):
        tsv_content = (
            "variant_id\tevidence\taf\n"
            "VAR1\tPVS1, PS3\t0.00001\n"
            "VAR2\tBA1\t0.15\n"
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write(tsv_content)
            temp_path = f.name

        try:
            out = io.StringIO()
            old_stdout = sys.stdout
            sys.stdout = out
            try:
                exit_code = cli.main(["-i", temp_path, "--format", "json"])
                self.assertEqual(exit_code, 0)
            finally:
                sys.stdout = old_stdout

            data = json.loads(out.getvalue())
            self.assertEqual(len(data), 2)
            self.assertEqual(data[0]["final_classification"], "Pathogenic")
            self.assertEqual(data[1]["final_classification"], "Benign")
        finally:
            if os.path.exists(temp_path):
                os.remove(temp_path)

    def test_cli_splice_mode(self):
        out = io.StringIO()
        old_stdout = sys.stdout
        sys.stdout = out
        try:
            exit_code = cli.main(["--splice", "--variant-id", "SP-01", "--pos", "1", "--ref", "G", "--alt", "A", "--format", "json"])
            self.assertEqual(exit_code, 0)
        finally:
            sys.stdout = old_stdout

        data = json.loads(out.getvalue())
        self.assertEqual(data["variant_id"], "SP-01")
        self.assertEqual(data["evidence_level"], "strong")

    def test_cli_pgx_mode(self):
        out = io.StringIO()
        old_stdout = sys.stdout
        sys.stdout = out
        try:
            exit_code = cli.main(["--pgx", "--gene", "CYP2C19", "--variant-id", "c.681G>A", "--format", "json"])
            self.assertEqual(exit_code, 0)
        finally:
            sys.stdout = old_stdout

        data = json.loads(out.getvalue())
        self.assertEqual(data["gene"], "CYP2C19")
        self.assertEqual(data["cpic_level"], "A")


if __name__ == "__main__":
    unittest.main()
