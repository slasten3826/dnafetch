# translate.py — TRANSLATE layer
# Maps DNA nucleotides to ProcessLang operators
# Calculates all statistical metrics for .dpl generation
# Part of DNAfetch v3

import math
from dataclasses import dataclass, field
from typing import List, Dict, Tuple
from collections import Counter

from fetch import Gene, AA_PROPERTIES

# ============================================================
# CORE MAPPING: Nucleotide → ProcessLang Operator
# ============================================================
# A = FLOW     (FluidCore)          — adenine, energy currency (ATP)
# T = OBSERVE  (SelfObservation)    — thymine, complementary witness
# G = LOGIC    (LogicSimulator)     — guanine, triple bond = structure
# C = CHOOSE   (ChoiceMechanism)    — cytosine, collapse of possibilities

NUC_TO_OP = {
    "A": "FLOW",
    "T": "OBSERVE",
    "G": "LOGIC",
    "C": "CHOOSE"
}

OP_TO_NUC = {v: k for k, v in NUC_TO_OP.items()}

# Codon positions
POS_NAMES = ["primary", "modifier", "wobble"]


@dataclass
class PositionProfile:
    """Operator frequencies at a specific codon position."""
    position: str  # primary / modifier / wobble
    counts: Dict[str, int] = field(default_factory=lambda: {
        "FLOW": 0, "OBSERVE": 0, "LOGIC": 0, "CHOOSE": 0
    })
    total: int = 0

    def add(self, operator: str):
        self.counts[operator] = self.counts.get(operator, 0) + 1
        self.total += 1

    @property
    def frequencies(self) -> Dict[str, float]:
        if self.total == 0:
            return {k: 0.0 for k in self.counts}
        return {k: v / self.total for k, v in self.counts.items()}

    @property
    def dominant(self) -> str:
        return max(self.counts, key=self.counts.get)

    @property
    def dominant_pct(self) -> float:
        if self.total == 0:
            return 0.0
        return self.counts[self.dominant] / self.total

    @property
    def entropy(self) -> float:
        """Shannon entropy of operator distribution (0 = uniform, 1 = one dominates)."""
        freqs = self.frequencies
        h = 0.0
        for f in freqs.values():
            if f > 0:
                h -= f * math.log2(f)
        max_h = math.log2(4)  # max entropy for 4 symbols
        return 1.0 - (h / max_h) if max_h > 0 else 0.0  # inverted: 1 = biased


@dataclass
class ComposePair:
    """A COMPOSE pair: two adjacent codon positions or codons."""
    op1: str
    op2: str
    count: int = 0
    frequency: float = 0.0

    @property
    def label(self) -> str:
        return f"COMPOSE({self.op1},{self.op2})"


@dataclass
class TranslationResult:
    """Complete ProcessLang translation of a gene."""
    gene: Gene

    # Position profiles
    primary: PositionProfile = None
    modifier: PositionProfile = None
    wobble: PositionProfile = None

    # Overall operator balance
    overall: PositionProfile = None

    # COMPOSE pairs (primary→modifier within codons)
    compose_pairs: List[ComposePair] = field(default_factory=list)

    # Derived metrics
    gc_content: float = 0.0
    wobble_gc: float = 0.0
    complementary_balance: float = 0.0  # how close A≈T, G≈C
    top2_gap: float = 0.0              # gap between top 2 operators overall
    codon_diversity: float = 0.0       # unique codons / 64
    hydrophobic_ratio: float = 0.0
    charge_profile: Dict[str, float] = field(default_factory=dict)
    aa_distribution: Dict[str, float] = field(default_factory=dict)
    aa_std: float = 0.0


def translate_gene(gene: Gene) -> TranslationResult:
    """Translate a Gene into ProcessLang operators and compute all metrics."""

    result = TranslationResult(gene=gene)

    # Initialize position profiles
    result.primary = PositionProfile(position="primary")
    result.modifier = PositionProfile(position="modifier")
    result.wobble = PositionProfile(position="wobble")
    result.overall = PositionProfile(position="overall")

    profiles = [result.primary, result.modifier, result.wobble]

    # Count operators by position
    for codon in gene.codons:
        for i, nuc in enumerate(codon):
            if nuc in NUC_TO_OP:
                op = NUC_TO_OP[nuc]
                profiles[i].add(op)
                result.overall.add(op)

    # COMPOSE pairs (primary → modifier)
    compose_counts = Counter()
    for codon in gene.codons:
        if len(codon) == 3:
            op1 = NUC_TO_OP.get(codon[0], "?")
            op2 = NUC_TO_OP.get(codon[1], "?")
            if op1 != "?" and op2 != "?":
                compose_counts[(op1, op2)] += 1

    total_compose = sum(compose_counts.values())
    result.compose_pairs = sorted([
        ComposePair(
            op1=k[0], op2=k[1],
            count=v,
            frequency=v / total_compose if total_compose > 0 else 0.0
        )
        for k, v in compose_counts.items()
    ], key=lambda x: -x.count)

    # GC content
    cds = gene.cds
    gc_count = sum(1 for n in cds if n in "GC")
    result.gc_content = gc_count / len(cds) if cds else 0.0

    # Wobble GC
    wobble_nucs = [codon[2] for codon in gene.codons if len(codon) == 3]
    wobble_gc = sum(1 for n in wobble_nucs if n in "GC")
    result.wobble_gc = wobble_gc / len(wobble_nucs) if wobble_nucs else 0.0

    # Complementary balance: |A-T| + |G-C| normalized
    of = result.overall.frequencies
    at_diff = abs(of.get("FLOW", 0) - of.get("OBSERVE", 0))
    gc_diff = abs(of.get("LOGIC", 0) - of.get("CHOOSE", 0))
    result.complementary_balance = 1.0 - (at_diff + gc_diff)

    # Top 2 gap
    sorted_ops = sorted(of.values(), reverse=True)
    result.top2_gap = (sorted_ops[0] - sorted_ops[1]) if len(sorted_ops) >= 2 else 0.0

    # Codon diversity
    unique_codons = len(set(gene.codons))
    result.codon_diversity = unique_codons / 64.0

    # Amino acid analysis
    aa_counts = Counter(aa for aa in gene.amino_acids if aa != "STOP")
    total_aa = sum(aa_counts.values())

    if total_aa > 0:
        result.aa_distribution = {
            aa: count / total_aa for aa, count in aa_counts.most_common()
        }

        # Hydrophobic ratio
        hydrophobic = sum(
            count for aa, count in aa_counts.items()
            if AA_PROPERTIES.get(aa, {}).get("hydrophobic", False)
        )
        result.hydrophobic_ratio = hydrophobic / total_aa

        # Charge profile
        charges = {"positive": 0, "negative": 0, "neutral": 0}
        for aa, count in aa_counts.items():
            ch = AA_PROPERTIES.get(aa, {}).get("charge", "neutral")
            charges[ch] += count
        result.charge_profile = {k: v / total_aa for k, v in charges.items()}

        # AA standard deviation (evenness of distribution)
        expected = 1.0 / 20.0
        variance = sum((f - expected) ** 2 for f in result.aa_distribution.values())
        # account for missing AAs
        missing = 20 - len(result.aa_distribution)
        variance += missing * (expected ** 2)
        result.aa_std = math.sqrt(variance / 20)

    return result
