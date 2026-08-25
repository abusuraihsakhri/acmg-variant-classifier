"""Command-line interface for the ACMG Variant Classifier."""

import argparse
import csv
import re
import sys

from . import report as report_mod
from .classifier import classify_variant
from .frequency import DEFAULT_BA1_THRESHOLD, DEFAULT_PM2_THRESHOLD


def _split_codes(text):
    if not text:
        return []
    return [c for c in re.split(r"[,;|\s]+", text.strip()) if c]


def build_parser():
    parser = argparse.ArgumentParser(
        prog="acmg-classify",
        description=(
            "Classify germline sequence variants as Pathogenic, Likely "
            "Pathogenic, Uncertain Significance, Likely Benign, or Benign "
            "using the ACMG/AMP 2015 combining rules and the Tavtigian "
            "2020 Bayesian point-scoring refinement."
        ),
    )

    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--evidence", "-e", nargs="+", metavar="CODE",
        help="Evidence codes for a single variant, e.g. -e PVS1 PS3 PM2 "
             "(codes may also be comma-separated).",
    )
    mode.add_argument(
        "--input", "-i", metavar="TSV_FILE",
        help="Batch mode: TSV file with columns 'variant_id', 'evidence' "
             "(codes separated by , ; or |), and optional 'af'.",
    )
    mode.add_argument(
        "--vcf", metavar="VCF_FILE",
        help="Batch mode: VCF file where the INFO column carries an "
             "'ACMG=' key (pipe/comma separated evidence codes) and "
             "optionally an 'AF=' key (population allele frequency).",
    )

    parser.add_argument(
        "--variant-id", default="variant",
        help="Variant identifier to use in the report (single-variant mode only).",
    )
    parser.add_argument(
        "--af", type=float, default=None,
        help="Population allele frequency for the variant (single-variant mode only).",
    )
    parser.add_argument(
        "--ba1-threshold", type=float, default=DEFAULT_BA1_THRESHOLD,
        help=f"Allele frequency above which BA1 (stand-alone benign) applies "
             f"(default: {DEFAULT_BA1_THRESHOLD}).",
    )
    parser.add_argument(
        "--pm2-threshold", type=float, default=DEFAULT_PM2_THRESHOLD,
        help=f"Allele frequency at/below which PM2 (absent/rare) applies "
             f"(default: {DEFAULT_PM2_THRESHOLD}).",
    )
    parser.add_argument(
        "--no-auto-frequency", action="store_true",
        help="Do not auto-derive BA1/PM2 from the population allele "
             "frequency; only report contradiction warnings.",
    )
    parser.add_argument(
        "--format", choices=["text", "json"], default="text",
        help="Output format (default: text).",
    )
    parser.add_argument(
        "--output", "-o", metavar="FILE",
        help="Write the report to FILE instead of stdout.",
    )
    return parser


def _classify_one(codes, variant_id, af, args):
    return classify_variant(
        codes,
        variant_id=variant_id,
        population_af=af,
        ba1_threshold=args.ba1_threshold,
        pm2_threshold=args.pm2_threshold,
        auto_frequency=not args.no_auto_frequency,
    )


def _run_single(args):
    codes = []
    for item in args.evidence:
        codes.extend(_split_codes(item))
    return [_classify_one(codes, args.variant_id, args.af, args)]


def _run_tsv_batch(args):
    reports = []
    with open(args.input, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None or "evidence" not in reader.fieldnames:
            raise ValueError(
                "TSV input must have a header row with at least a "
                "'evidence' column (and optionally 'variant_id', 'af')."
            )
        for row_num, row in enumerate(reader, start=2):
            variant_id = (row.get("variant_id") or f"row{row_num}").strip()
            codes = _split_codes(row.get("evidence", ""))
            af_raw = (row.get("af") or "").strip()
            af = float(af_raw) if af_raw else None
            reports.append(_classify_one(codes, variant_id, af, args))
    return reports


_VCF_INFO_ACMG_RE = re.compile(r"(?:^|;)ACMG=([^;]+)")
_VCF_INFO_AF_RE = re.compile(r"(?:^|;)AF=([^;]+)")


def _run_vcf_batch(args):
    reports = []
    with open(args.vcf, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, vid, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
            info = fields[7]

            acmg_match = _VCF_INFO_ACMG_RE.search(info)
            if not acmg_match:
                continue
            codes = _split_codes(acmg_match.group(1).replace("|", ","))

            af = None
            af_match = _VCF_INFO_AF_RE.search(info)
            if af_match:
                try:
                    af = float(af_match.group(1).split(",")[0])
                except ValueError:
                    af = None

            variant_id = vid if vid and vid != "." else f"{chrom}:{pos}:{ref}>{alt}"
            reports.append(_classify_one(codes, variant_id, af, args))
    return reports


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        if args.evidence is not None:
            reports = _run_single(args)
        elif args.input is not None:
            reports = _run_tsv_batch(args)
        else:
            reports = _run_vcf_batch(args)
    except (ValueError, FileNotFoundError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    if not reports:
        print("No variants were classified (no matching evidence found).", file=sys.stderr)
        return 1

    payload = reports if len(reports) > 1 else reports[0]
    if args.format == "json":
        output = report_mod.reports_to_json(payload)
    else:
        output = report_mod.reports_to_text(payload)

    if args.output:
        with open(args.output, "w", encoding="utf-8") as fh:
            fh.write(output + "\n")
    else:
        print(output)

    return 0


if __name__ == "__main__":
    sys.exit(main())
