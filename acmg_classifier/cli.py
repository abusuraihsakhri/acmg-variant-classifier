#!/usr/bin/env python3
"""Command Line Interface for ACMG/AMP Variant Classifier & Functional Genomics.

Supports interactive mode, single-variant evaluation, TSV/CSV/VCF batch processing,
splice impact prediction, and pharmacogenomics (CPIC) annotation.
"""

import argparse
import csv
import json
import math
import os
import re
import sys
from pathlib import Path
from typing import List, Optional, Union

from . import report as report_mod
from .classifier import classify_variant, ClassificationReport
from .frequency import DEFAULT_BA1_THRESHOLD, DEFAULT_PM2_THRESHOLD
from .security import safe_resolve_path, sanitize_csv_cell
from variant_analysis import SpliceImpactPredictor, PharmacoGenomicsAnnotator, FunctionalDomainMapper


def _split_codes(text: str) -> List[str]:
    """Split comma, semicolon, pipe, or whitespace-separated evidence codes."""
    if not text:
        return []
    return [c.strip() for c in re.split(r"[,;|\s]+", text.strip()) if c.strip()]


def _validate_input_path(filepath: Union[str, Path], base_dir: Optional[Union[str, Path]] = None) -> Path:
    """Validate that an input file path exists and is a regular file safely."""
    path = safe_resolve_path(filepath, base_dir=base_dir)
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {filepath}")
    if not path.is_file():
        raise ValueError(f"Input path is not a regular file: {filepath}")
    return path


def _validate_output_path(filepath: Union[str, Path], base_dir: Optional[Union[str, Path]] = None) -> Path:
    """Validate that an output file path is safe to write."""
    path = safe_resolve_path(filepath, base_dir=base_dir)
    if path.exists() and path.is_dir():
        raise ValueError(f"Output path is an existing directory: {filepath}")
    parent = path.parent
    if not parent.exists():
        parent.mkdir(parents=True, exist_ok=True)
    return path


def _probability_float(value: str) -> float:
    """Argparse type validator for probability values in [0, 1].

    Rejects non-numeric tokens, NaN, and Infinity.
    """
    if isinstance(value, bool):
        raise argparse.ArgumentTypeError(f"Invalid float value: {value!r}")
    try:
        fval = float(value)
    except (ValueError, TypeError):
        raise argparse.ArgumentTypeError(f"Invalid float value: {value!r}")
    if math.isnan(fval) or math.isinf(fval):
        raise argparse.ArgumentTypeError(f"Value cannot be NaN or infinite: {value!r}")
    if not 0.0 <= fval <= 1.0:
        raise argparse.ArgumentTypeError(
            f"Value {fval} is outside the valid probability range [0, 1]."
        )
    return fval


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="acmg-classifier",
        description="Production ACMG/AMP 2015 & Tavtigian 2020 Bayesian Variant Classifier & Genomic Annotator",
    )

    # Primary operation modes
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--evidence", "-e", nargs="+", metavar="CODE",
        help="ACMG evidence codes for a single variant (e.g., -e PVS1 PS3 PM2 or -e 'PVS1, PS3, PM2').",
    )
    mode.add_argument(
        "--input", "-i", metavar="FILE",
        help="Batch mode: CSV or TSV file containing 'variant_id', 'evidence', and optional 'af'.",
    )
    mode.add_argument(
        "--vcf", metavar="VCF_FILE",
        help="Batch mode: VCF file with INFO field tags 'ACMG=' and optional 'AF='.",
    )
    mode.add_argument(
        "--interactive", "-I", action="store_true",
        help="Launch interactive terminal wizard for step-by-step variant evaluation.",
    )
    mode.add_argument(
        "--splice", action="store_true",
        help="Predict in-silico splice site impact using position-weight matrices.",
    )
    mode.add_argument(
        "--pgx", action="store_true",
        help="Annotate CPIC pharmacogenomics level and drug guidelines for a gene/variant.",
    )
    mode.add_argument(
        "--domain", action="store_true",
        help="Map protein position to functional domains and hotspot residues.",
    )

    # Variant / Gene parameters
    parser.add_argument("--variant-id", default="VAR-001", help="Variant identifier (e.g. 'BRCA1:c.5266dupC').")
    parser.add_argument("--gene", default="BRCA1", help="Gene symbol for PGx or domain mapping.")
    parser.add_argument("--af", type=_probability_float, default=None, help="Population allele frequency (e.g. 0.00004).")
    parser.add_argument("--pos", type=int, default=1, help="Genomic or protein position for splice/domain mapping.")
    parser.add_argument("--ref", default="G", help="Reference base for splice prediction.")
    parser.add_argument("--alt", default="A", help="Alternate base for splice prediction.")
    parser.add_argument("--distance", type=int, default=0, help="Distance from splice junction in nucleotides.")
    parser.add_argument("--region-type", choices=["intron", "exon"], default="intron", help="Region type for splice prediction.")
    parser.add_argument("--variant-type", default="missense", help="Variant type (missense, nonsense, frameshift, splice).")

    # Algorithm thresholds
    parser.add_argument("--ba1-threshold", type=_probability_float, default=DEFAULT_BA1_THRESHOLD,
                        help=f"Allele frequency above which BA1 applies (default: {DEFAULT_BA1_THRESHOLD}).")
    parser.add_argument("--pm2-threshold", type=_probability_float, default=DEFAULT_PM2_THRESHOLD,
                        help=f"Allele frequency at/below which PM2 applies (default: {DEFAULT_PM2_THRESHOLD}).")
    parser.add_argument("--no-auto-frequency", action="store_true",
                        help="Disable auto-derivation of BA1/PM2 from allele frequency.")

    # Output formatting
    parser.add_argument("--format", choices=["text", "json", "csv"], default="text", help="Output format (default: text).")
    parser.add_argument("--json", action="store_true", help="Output result as formatted JSON (shorthand for --format json).")
    parser.add_argument("--output", "-o", metavar="FILE", help="Write report to specified output file.")

    return parser


def run_interactive(args) -> List[dict]:
    print("=" * 70)
    print("  ACMG/AMP 2015 & ClinGen Bayesian Variant Classification Wizard")
    print("=" * 70)

    try:
        variant_id = input("Enter Variant Identifier [default: VAR-001]: ").strip() or "VAR-001"
        evidence_input = input("Enter ACMG Evidence Codes (space or comma-separated, e.g. PVS1 PS3 PM2): ").strip()
        codes = _split_codes(evidence_input)

        af_input = input("Enter Population Allele Frequency (optional, e.g. 0.00002, leave blank if none): ").strip()
        af = float(af_input) if af_input else None

        gene = input("Enter Gene Symbol (optional, e.g. BRCA1, leave blank to skip): ").strip()
        pos_input = input("Enter Protein Position (optional, e.g. 64, leave blank to skip): ").strip()
        pos = int(pos_input) if pos_input else None

    except (KeyboardInterrupt, EOFError):
        print("\nInteractive mode aborted.")
        sys.exit(0)

    report = classify_variant(
        codes,
        variant_id=variant_id,
        population_af=af,
        ba1_threshold=args.ba1_threshold,
        pm2_threshold=args.pm2_threshold,
        auto_frequency=not args.no_auto_frequency,
    )

    print("\n" + "=" * 70)
    print(f"  CLASSIFICATION REPORT FOR: {variant_id}")
    print("=" * 70)
    print(report_mod.report_to_text(report))

    if gene and pos:
        mapper = FunctionalDomainMapper()
        hits = mapper.map_variant(gene, variant_id, pos)
        if hits:
            print("\n  FUNCTIONAL DOMAIN MAPPING:")
            for h in hits:
                print(f"  - Domain: {h.domain_name} ({h.residue_range[0]}-{h.residue_range[1]}) | In Domain: {h.in_domain} | Hot Spot: {h.hot_spot}")

    if gene:
        pgx = PharmacoGenomicsAnnotator()
        annot = pgx.annotate(gene, variant_id, args.variant_type)
        if annot.cpic_level != "No Data":
            print("\n  PHARMACOGENOMICS (CPIC):")
            print(f"  - Level: {annot.cpic_level} | Phenotype: {annot.phenotype}")
            print(f"  - Associated Drugs: {', '.join(annot.drug_names)}")
            print(f"  - Recommendation: {annot.recommendation}")
    print("=" * 70)

    return [report_mod.report_to_dict(report)]


def run_batch_table(args) -> List[ClassificationReport]:
    reports = []
    filepath = _validate_input_path(args.input)
    delimiter = "\t" if filepath.suffix.lower() in [".tsv", ".tab"] else ","

    with open(filepath, mode="r", encoding="utf-8-sig") as f:
        sample = f.read(2048)
        f.seek(0)
        if "\t" in sample and "," not in sample:
            delimiter = "\t"
        reader = csv.DictReader(f, delimiter=delimiter)

        for row_idx, row in enumerate(reader, start=1):
            vid = (row.get("variant_id") or row.get("variant") or f"row_{row_idx}").strip()
            ev_str = row.get("evidence") or row.get("codes") or ""
            codes = _split_codes(ev_str)
            af_val = row.get("af") or row.get("allele_frequency") or ""

            try:
                af = None
                if af_val.strip():
                    parsed_af = float(af_val.strip())
                    if math.isnan(parsed_af) or math.isinf(parsed_af) or not (0.0 <= parsed_af <= 1.0):
                        raise ValueError(f"Allele frequency {af_val!r} outside valid [0, 1] range")
                    af = parsed_af

                rep = classify_variant(
                    codes,
                    variant_id=vid,
                    population_af=af,
                    ba1_threshold=args.ba1_threshold,
                    pm2_threshold=args.pm2_threshold,
                    auto_frequency=not args.no_auto_frequency,
                )
                reports.append(rep)
            except (ValueError, TypeError) as exc:
                reports.append(ClassificationReport(
                    variant_id=vid,
                    input_codes=codes,
                    unknown_codes=[],
                    effective_codes=[],
                    auto_added_codes=[],
                    evidence_counts={},
                    categorical_label="Uncertain Significance",
                    categorical_triggered_rules=[],
                    categorical_conflict=False,
                    bayesian_points=0,
                    bayesian_label="Uncertain Significance",
                    final_classification="Uncertain Significance",
                    final_source="Batch Row Error",
                    warnings=[f"Row {row_idx} error: {exc}"],
                    rationale=f"Failed to classify variant due to row error: {exc}",
                ))
    return reports


def run_vcf(args) -> List[ClassificationReport]:
    reports = []
    vcf_acmg_re = re.compile(r"(?:^|;)ACMG=([^;]+)")
    vcf_af_re = re.compile(r"(?:^|;)AF=([^;]+)")

    vcf_path = _validate_input_path(args.vcf)
    with open(vcf_path, "r", encoding="utf-8") as fh:
        for line_idx, line in enumerate(fh, start=1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, vid, ref, alt, _, _, info = parts[:8]
            acmg_m = vcf_acmg_re.search(info)
            if not acmg_m:
                continue
            codes = _split_codes(acmg_m.group(1).replace("|", ","))
            af_m = vcf_af_re.search(info)
            af = None
            if af_m:
                try:
                    parsed_af = float(af_m.group(1).split(",")[0])
                    if not (math.isnan(parsed_af) or math.isinf(parsed_af) or parsed_af < 0.0 or parsed_af > 1.0):
                        af = parsed_af
                except (ValueError, TypeError):
                    af = None

            variant_id = vid if vid and vid != "." else f"{chrom}:{pos}:{ref}>{alt}"

            try:
                rep = classify_variant(
                    codes,
                    variant_id=variant_id,
                    population_af=af,
                    ba1_threshold=args.ba1_threshold,
                    pm2_threshold=args.pm2_threshold,
                    auto_frequency=not args.no_auto_frequency,
                )
                reports.append(rep)
            except (ValueError, TypeError) as exc:
                reports.append(ClassificationReport(
                    variant_id=variant_id,
                    input_codes=codes,
                    unknown_codes=[],
                    effective_codes=[],
                    auto_added_codes=[],
                    evidence_counts={},
                    categorical_label="Uncertain Significance",
                    categorical_triggered_rules=[],
                    categorical_conflict=False,
                    bayesian_points=0,
                    bayesian_label="Uncertain Significance",
                    final_classification="Uncertain Significance",
                    final_source="Batch Row Error",
                    warnings=[f"VCF line {line_idx} error: {exc}"],
                    rationale=f"Failed to classify variant due to line error: {exc}",
                ))
    return reports


def _run_main_logic(args) -> int:
    """Execute the main CLI logic with validated arguments."""
    if args.splice:
        predictor = SpliceImpactPredictor()
        res = predictor.predict(
            variant_id=args.variant_id,
            position=args.pos,
            ref_base=args.ref,
            alt_base=args.alt,
            intron_exon=args.region_type,
            splice_distance=args.distance,
        )
        if args.format == "json":
            out = json.dumps(res.__dict__, indent=2, allow_nan=False)
        else:
            out = (
                f"Splice Impact Prediction for {res.variant_id}:\n"
                f"  Site Type:         {res.splice_site_type}\n"
                f"  Predicted Effect:  {res.predicted_effect}\n"
                f"  Delta Score:       {res.delta_score} (Conservation: {res.conservation_score})\n"
                f"  Evidence Level:    {res.evidence_level}\n"
            )
        if args.output:
            _validate_output_path(args.output).write_text(out, encoding="utf-8")
        else:
            print(out)
        return 0

    if args.pgx:
        annotator = PharmacoGenomicsAnnotator()
        annot = annotator.annotate(args.gene, args.variant_id, args.variant_type)
        if args.format == "json":
            out = json.dumps(annot.__dict__, indent=2, allow_nan=False)
        else:
            out = (
                f"CPIC Pharmacogenomic Annotation for {annot.gene} ({annot.variant_id}):\n"
                f"  CPIC Level:     {annot.cpic_level}\n"
                f"  Phenotype:      {annot.phenotype}\n"
                f"  Target Drugs:   {', '.join(annot.drug_names) or 'None'}\n"
                f"  Recommendation: {annot.recommendation}\n"
                f"  Guideline URL:  {annot.evidence_url}\n"
            )
        if args.output:
            _validate_output_path(args.output).write_text(out, encoding="utf-8")
        else:
            print(out)
        return 0

    if args.domain:
        mapper = FunctionalDomainMapper()
        hits = mapper.map_variant(args.gene, args.variant_id, args.pos)
        if args.format == "json":
            out = json.dumps([h.__dict__ for h in hits], indent=2, allow_nan=False)
        else:
            out = f"Domain Mapping for {args.gene} at position {args.pos}:\n"
            for h in hits:
                out += f"  - Domain: {h.domain_name} ({h.residue_range[0]}-{h.residue_range[1]}) | In Domain: {h.in_domain} | Hot Spot: {h.hot_spot}\n"
        if args.output:
            _validate_output_path(args.output).write_text(out, encoding="utf-8")
        else:
            print(out)
        return 0

    if args.evidence:
        codes = []
        for item in args.evidence:
            codes.extend(_split_codes(item))
        report = classify_variant(
            codes,
            variant_id=args.variant_id,
            population_af=args.af,
            ba1_threshold=args.ba1_threshold,
            pm2_threshold=args.pm2_threshold,
            auto_frequency=not args.no_auto_frequency,
        )
        reports = [report]

    elif args.input:
        reports = run_batch_table(args)

    elif args.vcf:
        reports = run_vcf(args)

    else:
        return 1

    if args.format == "json":
        output_str = report_mod.reports_to_json(reports)
    elif args.format == "csv":
        output_str = report_mod.reports_to_csv(reports)
    else:
        output_str = report_mod.reports_to_text(reports)

    if args.output:
        _validate_output_path(args.output).write_text(output_str, encoding="utf-8")
    else:
        print(output_str)

    return 0


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.json:
        args.format = "json"

    if args.interactive or (not args.evidence and not args.input and not args.vcf and not args.splice and not args.pgx and not args.domain):
        if len(sys.argv) == 1 or args.interactive:
            run_interactive(args)
            return 0
        parser.print_help()
        return 1

    try:
        return _run_main_logic(args)
    except (FileNotFoundError, ValueError, TypeError, PermissionError, OSError) as exc:
        if getattr(args, "json", False) or getattr(args, "format", "") == "json":
            print(json.dumps({"error": str(exc)}, indent=2))
        else:
            print(f"Error: {exc}", file=sys.stderr)
        return 1
    except Exception as exc:
        if getattr(args, "json", False) or getattr(args, "format", "") == "json":
            print(json.dumps({"error": f"Unexpected error: {exc}"}, indent=2))
        else:
            print(f"Unexpected error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
