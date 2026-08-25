"""Cross-checks a user-supplied population allele frequency against the
BA1 (stand-alone benign, variant too common to be disease-causing) and
PM2 (absent/rare in population controls) frequency thresholds.

When a population allele frequency is supplied, BA1 and PM2 are derived
automatically from it (unless auto-apply is disabled) and merged with the
manually supplied evidence codes. Any contradiction between the manually
asserted codes and the frequency data (e.g. PM2 claimed for a common
variant) is reported as a warning.
"""

DEFAULT_BA1_THRESHOLD = 0.05
DEFAULT_PM2_THRESHOLD = 0.0001


def check_frequency(codes, population_af, ba1_threshold=DEFAULT_BA1_THRESHOLD,
                     pm2_threshold=DEFAULT_PM2_THRESHOLD, auto_apply=True):
    """
    Args:
        codes: set/list of manually supplied evidence codes (validated).
        population_af: float allele frequency, or None if not supplied.
        ba1_threshold: AF above which BA1 (stand-alone benign) applies.
        pm2_threshold: AF at/below which PM2 (absent/rare) applies.
        auto_apply: if True, add BA1/PM2 to the effective code set based
            on the frequency data.

    Returns:
        (effective_codes: set[str], warnings: list[str], auto_added: list[str])
    """
    manual_codes = set(codes)
    effective = set(codes)
    warnings = []
    auto_added = []

    if population_af is None:
        return effective, warnings, auto_added

    if not 0.0 <= population_af <= 1.0:
        warnings.append(
            f"Population allele frequency {population_af!r} is outside the "
            f"valid range [0, 1]."
        )
        return effective, warnings, auto_added

    exceeds_ba1 = population_af > ba1_threshold
    within_pm2 = population_af <= pm2_threshold

    if auto_apply:
        if exceeds_ba1 and "BA1" not in effective:
            effective.add("BA1")
            auto_added.append("BA1")
        if within_pm2 and "PM2" not in effective:
            effective.add("PM2")
            auto_added.append("PM2")

    if exceeds_ba1 and "PM2" in manual_codes:
        warnings.append(
            f"Contradiction: PM2 (absent/rare in controls) was asserted, but "
            f"population AF {population_af:g} exceeds the BA1 threshold "
            f"({ba1_threshold:g}) -- this variant is common."
        )
    if within_pm2 and "BA1" in manual_codes:
        warnings.append(
            f"Contradiction: BA1 (benign, stand-alone) was asserted, but "
            f"population AF {population_af:g} is at or below the PM2 "
            f"threshold ({pm2_threshold:g}) -- this variant is rare."
        )
    if exceeds_ba1 and "BA1" not in manual_codes and not auto_apply:
        warnings.append(
            f"Population AF {population_af:g} exceeds the BA1 threshold "
            f"({ba1_threshold:g}) but BA1 was not asserted."
        )
    if within_pm2 and "PM2" not in manual_codes and not auto_apply:
        warnings.append(
            f"Population AF {population_af:g} is at or below the PM2 "
            f"threshold ({pm2_threshold:g}) but PM2 was not asserted."
        )

    return effective, warnings, auto_added
