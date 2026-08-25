"""Tavtigian et al. 2020 Bayesian point-scoring refinement.

Tavtigian SV, Harrison SM, Boucher KM, Biesecker LG. "Modeling the
ACMG/AMP variant classification guidelines as a Bayesian classification
framework." Genet Med. 2020;22(6):1041-1051.

The 2015 categorical evidence strengths are mapped onto points on a
log-odds-of-pathogenicity scale, using a base unit of 1 point (Supporting)
that doubles with each strength tier: Supporting=1, Moderate=2, Strong=4,
Very Strong=8 (and the mirror-image benign values, negated). Points are
summed across all applicable criteria and the total is mapped to one of
the five classification tiers via fixed thresholds. BA1 is excluded from
the sum and treated as a stand-alone, decisive benign criterion, matching
its treatment in the original 2015 guidelines.
"""

POINTS_PER_CATEGORY = {
    "PVS": 8,
    "PS": 4,
    "PM": 2,
    "PP": 1,
    "BS": -4,
    "BP": -1,
}

BAYESIAN_POINTS = POINTS_PER_CATEGORY

# Point thresholds from Tavtigian et al. 2020, Table 3.
PATHOGENIC_MIN = 10
LIKELY_PATHOGENIC_RANGE = (6, 9)
VUS_RANGE = (0, 5)
LIKELY_BENIGN_RANGE = (-6, -1)
BENIGN_MAX = -7


def compute_points(counts):
    """Sum Tavtigian points across PVS/PS/PM/PP/BS/BP counts. BA1 is
    excluded (handled as a stand-alone criterion elsewhere)."""
    return sum(POINTS_PER_CATEGORY[cat] * counts[cat] for cat in POINTS_PER_CATEGORY)


def classify_points(points):
    """Map a summed point total to (rank, label) using the Tavtigian 2020
    thresholds."""
    if points >= PATHOGENIC_MIN:
        return 2, "Pathogenic"
    if LIKELY_PATHOGENIC_RANGE[0] <= points <= LIKELY_PATHOGENIC_RANGE[1]:
        return 1, "Likely Pathogenic"
    if VUS_RANGE[0] <= points <= VUS_RANGE[1]:
        return 0, "Uncertain Significance"
    if LIKELY_BENIGN_RANGE[0] <= points <= LIKELY_BENIGN_RANGE[1]:
        return -1, "Likely Benign"
    return -2, "Benign"


def bayesian_classification(counts, ba1_present):
    """Return a dict describing the Bayesian point-based classification.

    If BA1 is present it is decisive (stand-alone benign) regardless of
    the point sum, mirroring the categorical treatment of BA1.
    """
    points = compute_points(counts)
    if ba1_present:
        rank, label = -2, "Benign"
    else:
        rank, label = classify_points(points)
    return {
        "points": points,
        "rank": rank,
        "label": label,
        "ba1_stand_alone": ba1_present,
    }


def posterior_probability(points: int, prior: float = 0.10, odds_base: float = 2.08) -> float:
    """Compute Tavtigian 2020 Bayesian posterior probability of pathogenicity.

    Args:
        points: Sum of Tavtigian evidence points.
        prior: Prior probability of pathogenicity (default 0.10).
        odds_base: Exponential base per point (default 2.08 for ln(OddsPath) unit).

    Returns:
        Posterior probability in [0.0, 1.0].
    """
    if prior <= 0.0:
        return 0.0
    if prior >= 1.0:
        return 1.0
    prior_odds = prior / (1.0 - prior)
    lr = odds_base ** float(points)
    post_odds = prior_odds * lr
    post_prob = post_odds / (1.0 + post_odds)
    return max(0.0, min(1.0, post_prob))
