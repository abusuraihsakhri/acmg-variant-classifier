import csv
import io
import json
from dataclasses import asdict

from .security import sanitize_csv_cell


def report_to_dict(report):
    return asdict(report)


def reports_to_json(reports, indent=2):
    if isinstance(reports, list):
        data = [report_to_dict(r) for r in reports]
    else:
        data = report_to_dict(reports)
    return json.dumps(data, indent=indent, allow_nan=False)


def reports_to_csv(reports) -> str:
    """Serialize one or more ClassificationReports to CSV format.

    Neutralizes formula injection triggers (CWE-1236) by passing all cell
    values through sanitize_csv_cell.
    """
    if not isinstance(reports, list):
        reports = [reports]

    output = io.StringIO()
    writer = csv.writer(output, lineterminator="\n")

    headers = [
        "variant_id",
        "final_classification",
        "final_source",
        "categorical_label",
        "bayesian_points",
        "bayesian_label",
        "input_codes",
        "effective_codes",
        "auto_added_codes",
        "warnings",
        "rationale",
    ]
    writer.writerow([sanitize_csv_cell(h) for h in headers])

    for r in reports:
        row = [
            sanitize_csv_cell(r.variant_id),
            sanitize_csv_cell(r.final_classification),
            sanitize_csv_cell(r.final_source),
            sanitize_csv_cell(r.categorical_label),
            sanitize_csv_cell(r.bayesian_points),
            sanitize_csv_cell(r.bayesian_label),
            sanitize_csv_cell("; ".join(r.input_codes)),
            sanitize_csv_cell("; ".join(r.effective_codes)),
            sanitize_csv_cell("; ".join(r.auto_added_codes)),
            sanitize_csv_cell(" | ".join(r.warnings)),
            sanitize_csv_cell(r.rationale),
        ]
        writer.writerow(row)

    return output.getvalue()


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
