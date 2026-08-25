"""Combines the categorical ACMG/AMP 2015 rule table with the Tavtigian
2020 Bayesian point refinement into a single final classification, with
frequency cross-checking and internal-contradiction detection.
"""

from dataclasses import dataclass, field
from typing import List, Optional

from . import evidence as evidence_mod
from . import rules
from . import bayesian
from . import frequency

RANK_LABELS = {
    2: "Pathogenic",
    1: "Likely Pathogenic",
    0: "Uncertain Significance",
    -1: "Likely Benign",
    -2: "Benign",
}


@dataclass
class ClassificationReport:
    variant_id: str
    input_codes: List[str]
    unknown_codes: List[str]
    effective_codes: List[str]
    auto_added_codes: List[str]
    evidence_counts: dict
    categorical_label: str
    categorical_triggered_rules: List[str]
    categorical_conflict: bool
    bayesian_points: int
    bayesian_label: str
    final_classification: str
    final_source: str
    warnings: List[str] = field(default_factory=list)
    rationale: str = ""


def classify_variant(codes, variant_id="variant",
                      population_af=None,
                      ba1_threshold=frequency.DEFAULT_BA1_THRESHOLD,
                      pm2_threshold=frequency.DEFAULT_PM2_THRESHOLD,
                      auto_frequency=True):
    """Classify a single variant given its evidence codes.

    Args:
        codes: iterable of raw evidence code strings (e.g. ["PS1", "pm2"]).
        variant_id: label for the variant, used only in the report.
        population_af: optional float population allele frequency.
        ba1_threshold: AF above which BA1 (stand-alone benign) applies.
        pm2_threshold: AF at/below which PM2 (absent/rare) applies.
        auto_frequency: if True, derive BA1/PM2 automatically from
            population_af and merge them with the manually supplied codes.

    Returns:
        ClassificationReport
    """
    normalized_codes, unknown_codes = evidence_mod.validate_codes(codes)
    if unknown_codes:
        raise ValueError(
            f"Unknown ACMG evidence code(s): {', '.join(unknown_codes)}. "
            f"Valid codes: {', '.join(sorted(evidence_mod.VALID_CODES))}"
        )

    effective_codes, freq_warnings, auto_added = frequency.check_frequency(
        normalized_codes,
        population_af,
        ba1_threshold=ba1_threshold,
        pm2_threshold=pm2_threshold,
        auto_apply=auto_frequency,
    )

    counts = evidence_mod.count_by_category(effective_codes)

    cat_result = rules.categorical_classification(counts)
    bayes_result = bayesian.bayesian_classification(counts, ba1_present=counts["BA"] >= 1)

    warnings = list(freq_warnings)

    if cat_result["conflict"]:
        warnings.append(
            "Conflicting evidence: both pathogenic-tier and benign-tier "
            "categorical rules fired. Manual review is strongly "
            "recommended before reporting a final classification."
        )

    cat_rank = cat_result["rank"]
    bayes_rank = bayes_result["rank"]

    disagreement = cat_rank != 0 and bayes_rank != 0 and (cat_rank > 0) != (bayes_rank > 0)

    if cat_result["conflict"]:
        final_rank = 0
        final_source = "categorical rule conflict"
    elif disagreement:
        final_rank = 0
        final_source = "categorical/Bayesian directional disagreement"
        warnings.append(
            f"The 2015 categorical rules indicate '{cat_result['label']}' "
            f"while the Tavtigian 2020 point system ({bayes_result['points']} "
            f"points) indicates '{bayes_result['label']}' -- these disagree "
            f"in direction, so the variant is reported as Uncertain "
            f"Significance pending manual review."
        )
    else:
        # Same direction (or one of them is VUS/0): take whichever tier is
        # stronger. This lets the Bayesian point system act as a genuine
        # refinement, capable of resolving evidence combinations that the
        # discrete 2015 table does not explicitly cover.
        if abs(bayes_rank) > abs(cat_rank):
            final_rank = bayes_rank
            final_source = "Tavtigian 2020 Bayesian point refinement"
        else:
            final_rank = cat_rank
            final_source = "ACMG/AMP 2015 categorical rule table"

    final_label = RANK_LABELS[final_rank]

    rationale = _build_rationale(
        variant_id, counts, cat_result, bayes_result,
        final_label, final_source, auto_added, warnings,
    )

    return ClassificationReport(
        variant_id=variant_id,
        input_codes=normalized_codes,
        unknown_codes=unknown_codes,
        effective_codes=sorted(effective_codes),
        auto_added_codes=auto_added,
        evidence_counts=counts,
        categorical_label=cat_result["label"],
        categorical_triggered_rules=cat_result["triggered_rules"],
        categorical_conflict=cat_result["conflict"],
        bayesian_points=bayes_result["points"],
        bayesian_label=bayes_result["label"],
        final_classification=final_label,
        final_source=final_source,
        warnings=warnings,
        rationale=rationale,
    )


def _build_rationale(variant_id, counts, cat_result, bayes_result,
                      final_label, final_source, auto_added, warnings):
    lines = []
    summary = ", ".join(
        f"{cat}={counts[cat]}" for cat in ("PVS", "PS", "PM", "PP", "BA", "BS", "BP")
        if counts[cat]
    ) or "none"
    lines.append(f"Evidence tallied: {summary}.")
    if auto_added:
        lines.append(
            f"Auto-derived from population allele frequency: {', '.join(auto_added)}."
        )
    if cat_result["triggered_rules"]:
        lines.append(
            "Categorical rules (Richards et al. 2015) triggered: "
            + "; ".join(cat_result["triggered_rules"]) + "."
        )
    else:
        lines.append("No categorical (Richards et al. 2015) rule combination was met.")
    lines.append(
        f"Tavtigian 2020 Bayesian point total: {bayes_result['points']} "
        f"points -> {bayes_result['label']}."
    )
    lines.append(f"Final classification: {final_label} (determined by {final_source}).")
    if warnings:
        lines.append("Warnings: " + " | ".join(warnings))
    return " ".join(lines)
