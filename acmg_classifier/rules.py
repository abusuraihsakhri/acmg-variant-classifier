"""The ACMG/AMP 2015 categorical rule-combination table (Richards et al.,
Genet Med 17(5):405-424, 2015, Table 5).

Each evaluate_* function takes evidence counts per strength category and
returns (fired: bool, triggered_rule_descriptions: list[str]).
"""


def evaluate_pathogenic(n_pvs, n_ps, n_pm, n_pp):
    triggered = []
    if n_pvs >= 1 and n_ps >= 1:
        triggered.append("Pathogenic (i-a): 1 PVS1 AND >=1 PS")
    if n_pvs >= 1 and n_pm >= 2:
        triggered.append("Pathogenic (i-b): 1 PVS1 AND >=2 PM")
    if n_pvs >= 1 and n_pm >= 1 and n_pp >= 1:
        triggered.append("Pathogenic (i-c): 1 PVS1 AND 1 PM AND 1 PP")
    if n_pvs >= 1 and n_pp >= 2:
        triggered.append("Pathogenic (i-d): 1 PVS1 AND >=2 PP")
    if n_ps >= 2:
        triggered.append("Pathogenic (ii): >=2 PS")
    if n_ps >= 1 and n_pm >= 3:
        triggered.append("Pathogenic (iii-a): 1 PS AND >=3 PM")
    if n_ps >= 1 and n_pm >= 2 and n_pp >= 2:
        triggered.append("Pathogenic (iii-b): 1 PS AND 2 PM AND >=2 PP")
    if n_ps >= 1 and n_pm >= 1 and n_pp >= 4:
        triggered.append("Pathogenic (iii-c): 1 PS AND 1 PM AND >=4 PP")
    return bool(triggered), triggered


def evaluate_likely_pathogenic(n_pvs, n_ps, n_pm, n_pp):
    triggered = []
    if n_pvs >= 1 and n_pm >= 1:
        triggered.append("Likely Pathogenic (1): 1 PVS1 AND 1 PM")
    if n_ps >= 1 and 1 <= n_pm <= 2:
        triggered.append("Likely Pathogenic (2): 1 PS AND 1-2 PM")
    if n_ps >= 1 and n_pp >= 2:
        triggered.append("Likely Pathogenic (3): 1 PS AND >=2 PP")
    if n_pm >= 3:
        triggered.append("Likely Pathogenic (4): >=3 PM")
    if n_pm >= 2 and n_pp >= 2:
        triggered.append("Likely Pathogenic (5): 2 PM AND >=2 PP")
    if n_pm >= 1 and n_pp >= 4:
        triggered.append("Likely Pathogenic (6): 1 PM AND >=4 PP")
    return bool(triggered), triggered


def evaluate_benign(n_ba, n_bs):
    triggered = []
    if n_ba >= 1:
        triggered.append("Benign (i): 1 BA1 (stand-alone)")
    if n_bs >= 2:
        triggered.append("Benign (ii): >=2 BS")
    return bool(triggered), triggered


def evaluate_likely_benign(n_bs, n_bp):
    triggered = []
    if n_bs >= 1 and n_bp >= 1:
        triggered.append("Likely Benign (1): 1 BS AND 1 BP")
    if n_bp >= 2:
        triggered.append("Likely Benign (2): >=2 BP")
    return bool(triggered), triggered


def categorical_classification(counts):
    """Apply the full Richards et al. 2015 rule table to a dict of
    evidence counts (keys PVS, PS, PM, PP, BA, BS, BP).

    Returns a dict with:
      rank: int in {2, 1, 0, -1, -2} (2=Pathogenic .. -2=Benign, 0=VUS)
      label: str classification label
      triggered_rules: list[str] of every rule description that fired
      pathogenic_tier_hit: bool
      benign_tier_hit: bool
      conflict: bool (both a pathogenic-tier and benign-tier rule fired)
    """
    n_pvs, n_ps, n_pm, n_pp = counts["PVS"], counts["PS"], counts["PM"], counts["PP"]
    n_ba, n_bs, n_bp = counts["BA"], counts["BS"], counts["BP"]

    path_hit, path_rules = evaluate_pathogenic(n_pvs, n_ps, n_pm, n_pp)
    lp_hit, lp_rules = evaluate_likely_pathogenic(n_pvs, n_ps, n_pm, n_pp)
    benign_hit, benign_rules = evaluate_benign(n_ba, n_bs)
    lb_hit, lb_rules = evaluate_likely_benign(n_bs, n_bp)

    pathogenic_tier_hit = path_hit or lp_hit
    benign_tier_hit = benign_hit or lb_hit
    conflict = pathogenic_tier_hit and benign_tier_hit

    triggered_rules = path_rules + lp_rules + benign_rules + lb_rules

    if conflict:
        rank = 0
        label = "Uncertain Significance"
    elif path_hit:
        rank, label = 2, "Pathogenic"
    elif lp_hit:
        rank, label = 1, "Likely Pathogenic"
    elif benign_hit:
        rank, label = -2, "Benign"
    elif lb_hit:
        rank, label = -1, "Likely Benign"
    else:
        rank, label = 0, "Uncertain Significance"

    return {
        "rank": rank,
        "label": label,
        "triggered_rules": triggered_rules,
        "pathogenic_tier_hit": pathogenic_tier_hit,
        "benign_tier_hit": benign_tier_hit,
        "conflict": conflict,
    }
