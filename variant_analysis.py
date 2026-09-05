#!/usr/bin/env python3
"""
ACMG/AMP Variant Classifier: Splice Variant Impact Predictor & Pharmaco-Genomics Annotation.
Implements in-silico splice prediction, PGx CPIC level annotation, and phenotype-variant correlation.
"""
from dataclasses import dataclass
from typing import List, Tuple


@dataclass
class SplicePrediction:
    variant_id: str
    splice_site_type: str  # canonical, cryptic, exonic, intronic
    predicted_effect: str  # exon skipping, cryptic activation, no impact
    delta_score: float  # 0-1 severity score
    conservation_score: float
    evidence_level: str  # strong, moderate, supporting


@dataclass
class PGxAnnotation:
    variant_id: str
    gene: str
    cpic_level: str  # A, B, C, D, E, F, or "No Data"
    drug_names: List[str]
    phenotype: str
    recommendation: str
    evidence_url: str


@dataclass
class PhenotypeCorrelation:
    hpo_terms: List[str]
    zygosity: str  # heterozygous, homozygous, hemizygous
    inheritance: str
    penetrance: float
    expressivity: str


@dataclass
class FunctionalDomainHit:
    domain_name: str
    residue_range: Tuple[int, int]
    variant_position: int
    in_domain: bool
    hot_spot: bool


class SpliceImpactPredictor:
    """In-silico splice variant impact prediction using position-weight matrices."""

    CANONICAL_SPLICE_SITES = {
        "donor": {"gt": (0, 2), "ag": None},  # GT at +1,+2 of intron
        "acceptor": {"ag": (-2, -1), "ct": None},  # AG at -1,-2 of intron
    }

    CONSERVATION_POSITIONS = {
        "+1": "G", "+2": "T", "-1": "A", "-2": "G",
        "-3": "Y", "-4": "R", "+3": "R", "+4": "Y",
    }

    def predict(self, variant_id: str, position: int, ref_base: str,
                alt_base: str, intron_exon: str = "intron",
                splice_distance: int = 0) -> SplicePrediction:
        is_canonical = abs(splice_distance) <= 3 if intron_exon == "intron" else False
        delta_score = 0.0
        predicted_effect = "No splice impact predicted"
        site_type = "intronic"

        if intron_exon == "intron":
            if abs(splice_distance) <= 2:
                site_type = "canonical_donor" if splice_distance < 0 else "canonical_acceptor"
                delta_score = 0.9 if ref_base in ("G", "A", "T") else 0.7
                predicted_effect = "Loss of canonical splice site"
            elif abs(splice_distance) <= 7:
                site_type = "extended_splice_region"
                delta_score = 0.5
                predicted_effect = "Potential splice region disruption"
            else:
                site_type = "deep_intronic"
                delta_score = 0.2
                predicted_effect = "Possible cryptic splice site activation"
        elif intron_exon == "exon":
            site_type = "exonic_splice"
            delta_score = 0.3
            predicted_effect = "Possible exonic splice enhancer disruption"

        if alt_base in ("A", "T") and ref_base in ("C", "T"):
            delta_score = min(1.0, delta_score + 0.2)

        if delta_score >= 0.7:
            evidence = "strong"
        elif delta_score >= 0.4:
            evidence = "moderate"
        else:
            evidence = "supporting"

        conservation = self._assess_conservation(ref_base, position)

        return SplicePrediction(
            variant_id=variant_id,
            splice_site_type=site_type,
            predicted_effect=predicted_effect,
            delta_score=round(delta_score, 2),
            conservation_score=conservation,
            evidence_level=evidence,
        )

    def _assess_conservation(self, base: str, position: int) -> float:
        score = 0.5
        if position in (0, 1):
            score = 0.95
        elif position in (-1, -2):
            score = 0.90
        elif abs(position) <= 5:
            score = 0.7
        return score


class PharmacoGenomicsAnnotator:
    """ClinGen/CPIC pharmacogenomics annotation for variants."""

    CPIC_GENES = {
        "CYP2D6": {"drugs": ["codeine", "tamoxifen", "atomoxetine"], "level": "A",
                    "phenotype": "CYP2D6 poor/intermediate metabolizer"},
        "CYP2C19": {"drugs": ["clopidogrel", "omeprazole", "voriconazole"], "level": "A",
                     "phenotype": "CYP2C19 poor metabolizer"},
        "CYP2C9": {"drugs": ["warfarin", "phenytoin", "losartan"], "level": "A",
                    "phenotype": "CYP2C9 poor metabolizer"},
        "DPYD": {"drugs": ["5-fluorouracil", "capecitabine"], "level": "A",
                 "phenotype": "DPYD deficient metabolizer"},
        "TPMT": {"drugs": ["azathioprine", "6-mercaptopurine"], "level": "A",
                 "phenotype": "TPMT poor metabolizer"},
        "NUDT15": {"drugs": ["azathioprine", "6-mercaptopurine"], "level": "A",
                    "phenotype": "NUDT15 poor metabolizer"},
        "SLCO1B1": {"drugs": ["simvastatin", "atorvastatin"], "level": "A",
                     "phenotype": "SLCO1B1 poor function"},
        "VKORC1": {"drugs": ["warfarin"], "level": "A",
                   "phenotype": "VKORC1 low expression"},
        "HLA-B": {"drugs": ["carbamazepine", "abacavir", "allopurinol"], "level": "A",
                  "phenotype": "HLA-B allele-associated hypersensitivity"},
        "IFNL3": {"drugs": ["peginterferon alfa"], "level": "B",
                  "phenotype": "IFNL3 poor responder"},
    }

    def annotate(self, gene: str, variant_id: str, variant_type: str = "missense") -> PGxAnnotation:
        gene_upper = gene.upper()
        pgx_info = self.CPIC_GENES.get(gene_upper)

        if not pgx_info:
            return PGxAnnotation(
                variant_id=variant_id, gene=gene, cpic_level="No Data",
                drug_names=[], phenotype="No CPIC annotation available",
                recommendation="Consult PharmGKB for latest evidence",
                evidence_url="https://www.pharmgkb.org",
            )

        level = pgx_info["level"]
        if variant_type in ("frameshift", "nonsense", "splice"):
            recommendation = f"Likely affects {gene_upper} function. Review CPIC guidelines for {', '.join(pgx_info['drugs'][:2])}."
        elif variant_type == "missense":
            recommendation = f"May affect {gene_upper} function. Consider phenotyping if clinically relevant."
        else:
            recommendation = f"Assess functional impact of {variant_type} on {gene_upper}."

        return PGxAnnotation(
            variant_id=variant_id, gene=gene_upper, cpic_level=level,
            drug_names=pgx_info["drugs"], phenotype=pgx_info["phenotype"],
            recommendation=recommendation,
            evidence_url=f"https://www.pharmgkb.org/gene/{gene_upper}/clinicalAnnotation",
        )


class FunctionalDomainMapper:
    """Maps variant positions to protein functional domains for pathogenicity context."""

    DOMAIN_DB = {
        "BRCA1": [
            {"name": "RING domain", "start": 1, "end": 109, "hot_spot": True},
            {"name": "BRCT domain", "start": 1646, "end": 1863, "hot_spot": True},
            {"name": "Coiled-coil domain", "start": 1364, "end": 1437, "hot_spot": False},
        ],
        "TP53": [
            {"name": "DNA-binding domain", "start": 102, "end": 292, "hot_spot": True},
            {"name": "Tetramerization domain", "start": 325, "end": 356, "hot_spot": False},
            {"name": "Regulatory domain", "start": 1, "end": 60, "hot_spot": False},
        ],
        "CFTR": [
            {"name": "Nucleotide-binding domain 1", "start": 389, "end": 678, "hot_spot": True},
            {"name": "Nucleotide-binding domain 2", "start": 1218, "end": 1526, "hot_spot": True},
            {"name": "Regulatory (R) domain", "start": 679, "end": 836, "hot_spot": False},
        ],
        "BRAF": [
            {"name": "RAS-binding domain", "start": 156, "end": 294, "hot_spot": False},
            {"name": "Kinase domain", "start": 448, "end": 723, "hot_spot": True},
        ],
    }

    HOT_SPOT_POSITIONS = {
        "TP53": {175, 245, 248, 249, 273, 282},
        "KRAS": {12, 13, 61},
        "PIK3CA": {542, 545, 1047},
    }

    def map_variant(self, gene: str, variant_id: str, position: int) -> List[FunctionalDomainHit]:
        gene_upper = gene.upper()
        domains = self.DOMAIN_DB.get(gene_upper, [])
        hot_spots = self.HOT_SPOT_POSITIONS.get(gene_upper, set())

        hits = []
        for domain in domains:
            in_domain = domain["start"] <= position <= domain["end"]
            hot_spot = position in hot_spots if hot_spots else False
            hits.append(FunctionalDomainHit(
                domain_name=domain["name"],
                residue_range=(domain["start"], domain["end"]),
                variant_position=position,
                in_domain=in_domain,
                hot_spot=hot_spot,
            ))

        if not hits and hot_spots:
            hits.append(FunctionalDomainHit(
                domain_name="Unknown domain",
                residue_range=(0, 0),
                variant_position=position,
                in_domain=False,
                hot_spot=position in hot_spots,
            ))

        return hits
