"""Renders a ClassificationReport (or list of them) as JSON or plain text."""

import json
from dataclasses import asdict


def report_to_dict(report):
    return asdict(report)


def reports_to_json(reports, indent=2):
    if isinstance(reports, list):
        data = [report_to_dict(r) for r in reports]
    else:
        data = report_to_dict(reports)
    return json.dumps(data, indent=indent)


def report_to_text(report):
    lines = []
    lines.append(f"Variant: {report.variant_id}")
    lines.append(f"Final classification: {report.final_classification}")
    lines.append(f"  determined by: {report.final_source}")
    lines.append(f"Categorical (ACMG/AMP 2015): {report.categorical_label}")
    lines.append(f"Bayesian (Tavtigian 2020): {report.bayesian_label} "
                 f"({report.bayesian_points} points)")
    lines.append(f"Input evidence codes: {', '.join(report.input_codes) or 'none'}")
    if report.auto_added_codes:
        lines.append(f"Auto-added from allele frequency: {', '.join(report.auto_added_codes)}")
    counts_str = ", ".join(
        f"{cat}={count}" for cat, count in report.evidence_counts.items() if count
    ) or "none"
    lines.append(f"Evidence counts: {counts_str}")
    if report.categorical_triggered_rules:
        lines.append("Triggered categorical rules:")
        for rule in report.categorical_triggered_rules:
            lines.append(f"  - {rule}")
    if report.warnings:
        lines.append("Warnings:")
        for w in report.warnings:
            lines.append(f"  ! {w}")
    lines.append(f"Rationale: {report.rationale}")
    return "\n".join(lines)


def reports_to_text(reports):
    if isinstance(reports, list):
        return "\n\n".join(report_to_text(r) for r in reports)
    return report_to_text(reports)
