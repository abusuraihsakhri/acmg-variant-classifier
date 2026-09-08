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

    def test_cli_json_flag(self):
        out = io.StringIO()
        old_stdout = sys.stdout
        sys.stdout = out
        try:
            exit_code = cli.main(["--evidence", "PVS1", "PS3", "PM2", "--json", "--variant-id", "TEST-JSON"])
            self.assertEqual(exit_code, 0)
        finally:
            sys.stdout = old_stdout

        payload = json.loads(out.getvalue())
        self.assertEqual(payload[0]["final_classification"], "Pathogenic")
        self.assertEqual(payload[0]["variant_id"], "TEST-JSON")

    def test_benchmark_dataset_consistency(self):
        benchmark_file = PROJECT_ROOT / "benchmark_dataset.json"
        with open(benchmark_file, "r", encoding="utf-8") as f:
            bench = json.load(f)

        for item in bench.get("benchmark_variants", []):
            rep = classify_variant(
                item["evidence_codes"],
                variant_id=item["variant_id"],
                population_af=item.get("allele_frequency"),
            )
            self.assertEqual(
                rep.categorical_label,
                item["expected_categorical_classification"],
                f"Mismatch categorical for {item['variant_id']}",
            )
            self.assertEqual(
                rep.bayesian_points,
                item["expected_bayesian_points"],
                f"Mismatch bayesian points for {item['variant_id']}",
            )
            self.assertEqual(
                rep.bayesian_label,
                item["expected_bayesian_classification"],
                f"Mismatch bayesian label for {item['variant_id']}",
            )


class TestInputValidationAndSecurity(unittest.TestCase):
    """Tests for input validation, path safety, and error handling."""

    def test_validate_input_path_not_found(self):
        """Non-existent input file raises FileNotFoundError."""
        with self.assertRaises(FileNotFoundError):
            cli._validate_input_path("/nonexistent/path/file.csv")

    def test_validate_input_path_is_directory(self):
        """Directory path raises ValueError."""
        import tempfile
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaises(ValueError):
                cli._validate_input_path(tmpdir)

    def test_validate_output_path_is_directory(self):
        """Output path pointing to existing directory raises ValueError."""
        import tempfile
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaises(ValueError):
                cli._validate_output_path(tmpdir)

    def test_validate_output_path_creates_parent(self):
        """Output path creates parent directory if needed."""
        import tempfile
        import os
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = os.path.join(tmpdir, "subdir", "output.txt")
            result = cli._validate_output_path(output_path)
            self.assertTrue(os.path.isdir(os.path.dirname(result)))

    def test_probability_float_valid(self):
        """Valid probability values are accepted."""
        self.assertEqual(cli._probability_float("0.5"), 0.5)
        self.assertEqual(cli._probability_float("0"), 0.0)
        self.assertEqual(cli._probability_float("1"), 1.0)
        self.assertEqual(cli._probability_float("0.0001"), 0.0001)

    def test_probability_float_out_of_range(self):
        """Out-of-range probability values raise ArgumentTypeError."""
        import argparse
        with self.assertRaises(argparse.ArgumentTypeError):
            cli._probability_float("1.5")
        with self.assertRaises(argparse.ArgumentTypeError):
            cli._probability_float("-0.1")

    def test_probability_float_invalid(self):
        """Non-numeric values raise ArgumentTypeError."""
        import argparse
        with self.assertRaises(argparse.ArgumentTypeError):
            cli._probability_float("not_a_number")

    def test_cli_invalid_af_rejected(self):
        """CLI rejects out-of-range allele frequency."""
        import argparse
        with self.assertRaises(SystemExit):
            cli.main(["--evidence", "PVS1", "--af", "2.0"])

    def test_cli_invalid_ba1_threshold_rejected(self):
        """CLI rejects out-of-range BA1 threshold."""
        with self.assertRaises(SystemExit):
            cli.main(["--evidence", "PVS1", "--ba1-threshold", "-0.5"])

    def test_split_codes_strips_whitespace(self):
        """_split_codes strips whitespace from individual codes."""
        result = cli._split_codes(" PVS1 , PS3 , PM2 ")
        self.assertEqual(result, ["PVS1", "PS3", "PM2"])

    def test_cli_batch_missing_file(self):
        """CLI handles missing input file gracefully."""
        import io
        import sys
        old_stderr = sys.stderr
        sys.stderr = io.StringIO()
        try:
            exit_code = cli.main(["-i", "/nonexistent/file.csv"])
            self.assertEqual(exit_code, 1)
        finally:
            sys.stderr = old_stderr


class TestDefensiveHardeningAndSecurity(unittest.TestCase):
    """Rigorous tests for OWASP Top 10 hardening and defensive design patterns."""

    def test_csv_formula_injection_sanitization(self):
        from acmg_classifier.security import sanitize_csv_cell
        self.assertEqual(sanitize_csv_cell("=1+1"), "'=1+1")
        self.assertEqual(sanitize_csv_cell("+cmd|' /C calc'!A0"), "'+cmd|' /C calc'!A0")
        self.assertEqual(sanitize_csv_cell("@SUM(A1:A10)"), "'@SUM(A1:A10)")
        self.assertEqual(sanitize_csv_cell("\tDANGEROUS"), "'\tDANGEROUS")
        self.assertEqual(sanitize_csv_cell("\rDANGEROUS"), "'\rDANGEROUS")
        self.assertEqual(sanitize_csv_cell("-cmd|' /C calc'!A0"), "'-cmd|' /C calc'!A0")
        # Legitimate numbers preserved
        self.assertEqual(sanitize_csv_cell(-4.5), -4.5)
        self.assertEqual(sanitize_csv_cell("-4.5"), "-4.5")
        self.assertEqual(sanitize_csv_cell("-100"), "-100")
        self.assertEqual(sanitize_csv_cell(42), 42)
        self.assertEqual(sanitize_csv_cell("normal_text"), "normal_text")
        self.assertEqual(sanitize_csv_cell(None), "")

    def test_reports_to_csv_sanitizes_content(self):
        from acmg_classifier.report import reports_to_csv
        rep = classify_variant(["PVS1"], variant_id="=cmd|' /C calc'!A0")
        csv_out = reports_to_csv([rep])
        self.assertIn("'=cmd|' /C calc'!A0", csv_out)
        self.assertNotIn("\n=cmd|", csv_out)

    def test_path_traversal_windows_reserved(self):
        from acmg_classifier.security import safe_resolve_path
        for name in ["CON", "PRN", "AUX", "NUL", "COM1", "LPT1"]:
            with self.assertRaises(ValueError):
                safe_resolve_path(f"{name}.txt")
            with self.assertRaises(ValueError):
                safe_resolve_path(name)

    def test_path_traversal_null_bytes(self):
        from acmg_classifier.security import safe_resolve_path
        with self.assertRaises(ValueError):
            safe_resolve_path("test\0file.txt")

    def test_path_traversal_base_dir_confinement(self):
        from acmg_classifier.security import safe_resolve_path
        with tempfile.TemporaryDirectory() as base_dir:
            safe_inside = os.path.join(base_dir, "safe.txt")
            resolved = safe_resolve_path(safe_inside, base_dir=base_dir)
            self.assertEqual(resolved, Path(os.path.realpath(safe_inside)))

            escaping_path = os.path.join(base_dir, "..", "escape.txt")
            with self.assertRaises(PermissionError):
                safe_resolve_path(escaping_path, base_dir=base_dir)

    def test_numeric_validation_nan_inf_bool(self):
        from acmg_classifier.security import validate_numeric_range
        # Normal bounds succeed
        validate_numeric_range("val", 0.5, 0.0, 1.0)

        # Booleans rejected
        with self.assertRaises(TypeError):
            validate_numeric_range("val", True, 0.0, 1.0)
        with self.assertRaises(TypeError):
            validate_numeric_range("val", False, 0.0, 1.0)

        # NaN rejected
        with self.assertRaises(ValueError):
            validate_numeric_range("val", float("nan"), 0.0, 1.0)

        # Inf rejected
        with self.assertRaises(ValueError):
            validate_numeric_range("val", float("inf"), 0.0, 1.0)
        with self.assertRaises(ValueError):
            validate_numeric_range("val", float("-inf"), 0.0, 1.0)

        # Out of bounds rejected
        with self.assertRaises(ValueError):
            validate_numeric_range("val", 1.5, 0.0, 1.0)
        with self.assertRaises(ValueError):
            validate_numeric_range("val", -0.1, 0.0, 1.0)

    def test_posterior_probability_type_rejection(self):
        # Booleans and NaNs must be rejected
        with self.assertRaises(TypeError):
            bayesian.posterior_probability(points=True, prior=0.1)
        with self.assertRaises(TypeError):
            bayesian.posterior_probability(points=5, prior=True)
        with self.assertRaises(ValueError):
            bayesian.posterior_probability(points=float("nan"), prior=0.1)
        with self.assertRaises(ValueError):
            bayesian.posterior_probability(points=5, prior=float("nan"))
        with self.assertRaises(ValueError):
            bayesian.posterior_probability(points=5, prior=1.5)

    def test_frequency_check_type_rejection(self):
        # Boolean population_af should not bypass into BA1
        effective_codes, warnings, auto_added = frequency.check_frequency(
            ["PP1"], population_af=True
        )
        self.assertNotIn("BA1", effective_codes)
        self.assertTrue(any("must be numeric" in w for w in warnings))

        # NaN population_af should not bypass
        effective_codes, warnings, auto_added = frequency.check_frequency(
            ["PP1"], population_af=float("nan")
        )
        self.assertNotIn("BA1", effective_codes)
        self.assertTrue(any("outside the valid range" in w for w in warnings))

    def test_splice_predictor_type_rejection(self):
        predictor = SpliceImpactPredictor()
        with self.assertRaises(TypeError):
            predictor.predict("V1", position=True, ref_base="G", alt_base="A")
        with self.assertRaises(TypeError):
            predictor.predict("V1", position=1, ref_base="G", alt_base="A", splice_distance=True)
        with self.assertRaises(ValueError):
            predictor.predict("V1", position=1, ref_base="G", alt_base="A", intron_exon="invalid")

    def test_domain_mapper_type_rejection(self):
        mapper = FunctionalDomainMapper()
        with self.assertRaises(TypeError):
            mapper.map_variant("BRCA1", "V1", position=True)
        with self.assertRaises(TypeError):
            mapper.map_variant("BRCA1", "V1", position="first")

    def test_resilient_batch_csv_processing(self):
        csv_content = (
            "variant_id,evidence,af\n"
            "ROW1,PVS1,0.00001\n"
            "ROW2,UNKNOWN_CODE_ABC,0.00001\n"
            "ROW3,BA1,0.15\n"
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write(csv_content)
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
            self.assertEqual(len(data), 3)
            self.assertEqual(data[0]["variant_id"], "ROW1")
            self.assertEqual(data[0]["final_classification"], "Pathogenic")
            # Second row had error, captured without crashing
            self.assertEqual(data[1]["variant_id"], "ROW2")
            self.assertEqual(data[1]["final_source"], "Batch Row Error")
            self.assertTrue(any("UNKNOWN_CODE_ABC" in w for w in data[1]["warnings"]))
            # Third row succeeded
            self.assertEqual(data[2]["variant_id"], "ROW3")
            self.assertEqual(data[2]["final_classification"], "Benign")
        finally:
            if os.path.exists(temp_path):
                os.remove(temp_path)

    def test_cli_json_error_output(self):
        out = io.StringIO()
        old_stdout = sys.stdout
        sys.stdout = out
        try:
            exit_code = cli.main(["--evidence", "UNKNOWN_XYZ", "--json"])
            self.assertEqual(exit_code, 1)
        finally:
            sys.stdout = old_stdout

        err_json = json.loads(out.getvalue())
        self.assertIn("error", err_json)
        self.assertIn("UNKNOWN_XYZ", err_json["error"])

    def test_cli_csv_format_output(self):
        out = io.StringIO()
        old_stdout = sys.stdout
        sys.stdout = out
        try:
            exit_code = cli.main(["--evidence", "PVS1", "PS3", "--format", "csv", "--variant-id", "VAR-CSV"])
            self.assertEqual(exit_code, 0)
        finally:
            sys.stdout = old_stdout

        content = out.getvalue()
        self.assertIn("variant_id,final_classification", content)
        self.assertIn("VAR-CSV,Pathogenic", content)


if __name__ == "__main__":
    unittest.main()

