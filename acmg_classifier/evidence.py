"""The 26 standard ACMG/AMP 2015 evidence codes.

Each code is mapped to its evidence category (PVS, PS, PM, PP, BA, BS, BP),
its direction ("pathogenic" or "benign"), the Tavtigian et al. 2020 point
value used by the Bayesian refinement, and a short description.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class EvidenceDefinition:
    code: str
    category: str          # PVS, PS, PM, PP, BA, BS, BP
    direction: str          # "pathogenic" or "benign"
    points: int              # Tavtigian 2020 point value (signed)
    description: str


# Point values from Tavtigian et al. 2020 ("Modeling the ACMG/AMP variant
# classification guidelines as a Bayesian classification framework"),
# derived from a base exponential unit of 1 (Supporting) doubling per
# strength tier: Supporting=1, Moderate=2, Strong=4, Very Strong=8.
# BA1 is treated as a stand-alone, decisive criterion (as in the 2015
# guidelines) rather than as part of the additive point sum.
_CATEGORY_POINTS = {
    "PVS": 8,
    "PS": 4,
    "PM": 2,
    "PP": 1,
    "BA": None,  # stand-alone, not summed
    "BS": -4,
    "BP": -1,
}

_RAW_CODES = [
    ("PVS1", "PVS", "pathogenic",
     "Null variant (nonsense, frameshift, canonical +/-1 or 2 splice site, "
     "initiation codon, single or multi-exon deletion) in a gene where "
     "loss of function is a known mechanism of disease."),
    ("PS1", "PS", "pathogenic",
     "Same amino acid change as a previously established pathogenic "
     "variant regardless of nucleotide change."),
    ("PS2", "PS", "pathogenic",
     "De novo (with confirmed maternity and paternity) in a patient with "
     "the disease and no family history."),
    ("PS3", "PS", "pathogenic",
     "Well-established in vitro or in vivo functional studies show a "
     "damaging effect on protein function or splicing."),
    ("PS4", "PS", "pathogenic",
     "Prevalence of the variant in affected individuals is significantly "
     "increased compared with controls."),
    ("PM1", "PM", "pathogenic",
     "Located in a mutational hot spot or well-established functional "
     "domain without benign variation."),
    ("PM2", "PM", "pathogenic",
     "Absent, or at extremely low frequency if recessive, from controls "
     "in population frequency databases."),
    ("PM3", "PM", "pathogenic",
     "For recessive disorders, detected in trans with a pathogenic "
     "variant."),
    ("PM4", "PM", "pathogenic",
     "Protein length change due to an in-frame indel or stop-loss variant "
     "in a non-repeat region."),
    ("PM5", "PM", "pathogenic",
     "Novel missense change at a residue where a different missense "
     "change has been established as pathogenic."),
    ("PM6", "PM", "pathogenic",
     "Assumed de novo, but without confirmation of paternity and "
     "maternity."),
    ("PP1", "PP", "pathogenic",
     "Co-segregation with disease in multiple affected family members."),
    ("PP2", "PP", "pathogenic",
     "Missense variant in a gene with a low rate of benign missense "
     "variation and where missense variants are a common mechanism of "
     "disease."),
    ("PP3", "PP", "pathogenic",
     "Multiple lines of computational evidence support a deleterious "
     "effect on the gene or gene product."),
    ("PP4", "PP", "pathogenic",
     "Patient phenotype or family history is highly specific for a "
     "disease with a single genetic etiology."),
    ("PP5", "PP", "pathogenic",
     "Reputable source reports the variant as pathogenic, but the "
     "evidence is not available to the laboratory to perform an "
     "independent evaluation."),
    ("BA1", "BA", "benign",
     "Allele frequency is greater than 5% (or the supplied stand-alone "
     "threshold) in a population frequency database."),
    ("BS1", "BS", "benign",
     "Allele frequency is greater than expected for the disorder."),
    ("BS2", "BS", "benign",
     "Observed in healthy adult individuals for a disease with full "
     "penetrance expected at an early age."),
    ("BS3", "BS", "benign",
     "Well-established in vitro or in vivo functional studies show no "
     "damaging effect on protein function or splicing."),
    ("BS4", "BS", "benign",
     "Lack of segregation in affected members of a family."),
    ("BP1", "BP", "benign",
     "Missense variant in a gene for which primarily truncating variants "
     "are known to cause disease."),
    ("BP2", "BP", "benign",
     "Observed in trans with a pathogenic variant for a fully penetrant "
     "dominant disorder, or in cis with a pathogenic variant regardless "
     "of phase."),
    ("BP3", "BP", "benign",
     "In-frame indel in a repetitive region without a known function."),
    ("BP4", "BP", "benign",
     "Multiple lines of computational evidence suggest no impact on the "
     "gene or gene product."),
    ("BP5", "BP", "benign",
     "Variant found in a case with an alternate molecular basis for "
     "disease."),
    ("BP6", "BP", "benign",
     "Reputable source reports the variant as benign, but the evidence "
     "is not available to the laboratory to perform an independent "
     "evaluation."),
    ("BP7", "BP", "benign",
     "Synonymous variant for which splicing prediction algorithms "
     "predict no impact on the splice consensus sequence and the "
     "nucleotide is not highly conserved."),
]

EVIDENCE_CODES = {
    code: EvidenceDefinition(
        code=code,
        category=category,
        direction=direction,
        points=_CATEGORY_POINTS[category],
        description=description,
    )
    for code, category, direction, description in _RAW_CODES
}

VALID_CODES = frozenset(EVIDENCE_CODES.keys())

CATEGORIES_PATHOGENIC = ("PVS", "PS", "PM", "PP")
CATEGORIES_BENIGN = ("BA", "BS", "BP")


def validate_codes(codes):
    """Return (normalized_unique_codes, unknown_codes)."""
    normalized = []
    unknown = []
    seen = set()
    for raw in codes:
        code = raw.strip().upper()
        if not code:
            continue
        if code not in VALID_CODES:
            unknown.append(raw)
            continue
        if code not in seen:
            seen.add(code)
            normalized.append(code)
    return normalized, unknown


def count_by_category(codes):
    """Given an iterable of valid evidence codes, count occurrences per
    category (PVS, PS, PM, PP, BA, BS, BP)."""
    counts = {cat: 0 for cat in CATEGORIES_PATHOGENIC + CATEGORIES_BENIGN}
    for code in codes:
        counts[EVIDENCE_CODES[code].category] += 1
    return counts
