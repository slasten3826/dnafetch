# describe.py — DESCRIBE layer
# Converts TranslationResult into .dpl (DNA ProcessLang) format
# Uses ProcessLang operator names as section headers
# Part of DNAfetch v3

from translate import TranslationResult, NUC_TO_OP, OP_TO_NUC
from fetch import AA_PROPERTIES


def _observation_level(top2_gap: float) -> str:
    if top2_gap < 0.005:
        return "intense"
    elif top2_gap < 0.02:
        return "strong"
    elif top2_gap < 0.05:
        return "deep"
    elif top2_gap < 0.08:
        return "standard"
    elif top2_gap < 0.12:
        return "mild"
    else:
        return "basic"


def _cycle_stability(compose_pairs) -> str:
    if not compose_pairs:
        return "none"
    top_freq = compose_pairs[0].frequency
    if top_freq > 0.12:
        return "dominant"
    elif top_freq > 0.09:
        return "strong"
    elif top_freq > 0.07:
        return "moderate"
    else:
        return "distributed"


def _encoding_efficiency(wobble_gc: float) -> str:
    if wobble_gc > 0.75:
        return "engineered"
    elif wobble_gc > 0.60:
        return "optimized"
    elif wobble_gc > 0.45:
        return "moderate"
    elif wobble_gc > 0.30:
        return "natural"
    else:
        return "AT_rich"


def _dissolution_potential(wobble_entropy: float, comp_balance: float) -> str:
    score = wobble_entropy * 0.6 + (1.0 - comp_balance) * 0.4
    if score > 0.5:
        return "high"
    elif score > 0.3:
        return "moderate"
    elif score > 0.15:
        return "low"
    else:
        return "crystallized"


def _connection_type(hydro_ratio: float, charge: dict) -> str:
    positive = charge.get("positive", 0)
    negative = charge.get("negative", 0)
    if hydro_ratio > 0.55:
        return "membrane_embedded"
    elif positive > 0.15 or negative > 0.15:
        return "ionic_interaction"
    elif hydro_ratio < 0.35:
        return "soluble_mediator"
    else:
        return "balanced_interface"


def _manifestation_class(aa_dist: dict, protein_length: int) -> str:
    gly = aa_dist.get("Gly", 0)
    pro = aa_dist.get("Pro", 0)
    cys = aa_dist.get("Cys", 0)
    leu = aa_dist.get("Leu", 0)
    ser = aa_dist.get("Ser", 0)
    if gly > 0.15 and pro > 0.10:
        return "structural"
    elif cys > 0.04 and protein_length < 200:
        return "signaling_peptide"
    elif leu > 0.12:
        return "transmembrane"
    elif ser > 0.10:
        return "regulatory"
    else:
        return "enzymatic"


def generate_dpl(result: TranslationResult) -> str:
    gene = result.gene
    r = result

    wobble_class = _encoding_efficiency(r.wobble_gc)
    obs_level = _observation_level(r.top2_gap)
    cycle_stab = _cycle_stability(r.compose_pairs)
    dissolution = _dissolution_potential(r.wobble.entropy, r.complementary_balance)
    connection = _connection_type(r.hydrophobic_ratio, r.charge_profile)
    manifestation = _manifestation_class(r.aa_distribution, gene.protein_length)

    lines = []
    L = lines.append

    L(f"# DNAfetch v3 — DNA ProcessLang (.dpl)")
    L(f"# Gene: {gene.name} | {gene.organism}")
    L(f"# CDS: {gene.cds_length} bp ({gene.codon_count} codons)")
    L(f"")

    L(f"[IDENTITY]")
    L(f"gene = {gene.name}")
    L(f"organism = {gene.organism}")
    L(f"transcript = {gene.transcript_id}")
    L(f"cds_length = {gene.cds_length}")
    L(f"codons = {gene.codon_count}")
    L(f"protein_length = {gene.protein_length}")
    L(f"stop_codon = {gene.stop_codon}")
    L(f"")

    of = r.overall.frequencies
    flow_state = "flowing" if r.overall.entropy < 0.3 else "resistant"
    L(f"[FLOW]")
    L(f"engagement = {1.0 - r.overall.entropy:.3f}")
    L(f"resistance = {r.overall.entropy:.3f}")
    L(f"dominant = {r.overall.dominant}")
    L(f"dominant_pct = {r.overall.dominant_pct:.4f}")
    L(f"FLOW = {of.get('FLOW', 0):.4f}")
    L(f"OBSERVE = {of.get('OBSERVE', 0):.4f}")
    L(f"LOGIC = {of.get('LOGIC', 0):.4f}")
    L(f"CHOOSE = {of.get('CHOOSE', 0):.4f}")
    L(f"state = {flow_state}")
    L(f"")

    L(f"[CONNECT]")
    L(f"depth = {r.complementary_balance:.4f}")
    L(f"type = {connection}")
    L(f"hydrophobic_ratio = {r.hydrophobic_ratio:.4f}")
    for ch, val in r.charge_profile.items():
        L(f"charge_{ch} = {val:.4f}")
    L(f"")

    L(f"[DISSOLVE]")
    L(f"rigidity = {1.0 - r.wobble.entropy:.4f}")
    L(f"potential = {dissolution}")
    L(f"wobble_entropy = {r.wobble.entropy:.4f}")
    L(f"complementary_disruption = {1.0 - r.complementary_balance:.4f}")
    L(f"")

    L(f"[ENCODE]")
    L(f"type = {wobble_class}")
    L(f"loss = {1.0 - r.codon_diversity:.4f}")
    L(f"gc_content = {r.gc_content:.4f}")
    L(f"wobble_gc = {r.wobble_gc:.4f}")
    L(f"codon_diversity = {r.codon_diversity:.4f}")
    L(f"")

    wf = r.wobble.frequencies
    L(f"[CHOOSE]")
    L(f"method = collapse")
    L(f"dominant = {r.wobble.dominant}")
    L(f"dominant_pct = {r.wobble.dominant_pct:.4f}")
    L(f"wobble_FLOW = {wf.get('FLOW', 0):.4f}")
    L(f"wobble_OBSERVE = {wf.get('OBSERVE', 0):.4f}")
    L(f"wobble_LOGIC = {wf.get('LOGIC', 0):.4f}")
    L(f"wobble_CHOOSE = {wf.get('CHOOSE', 0):.4f}")
    L(f"entropy = {r.wobble.entropy:.4f}")
    L(f"")

    L(f"[OBSERVE]")
    L(f"distance = {r.top2_gap:.4f}")
    L(f"level = {obs_level}")
    L(f"complementary_balance = {r.complementary_balance:.4f}")
    L(f"")

    L(f"[CYCLE]")
    L(f"stability = {cycle_stab}")
    if r.compose_pairs:
        top = r.compose_pairs[0]
        L(f"top_compose = {top.label}")
        L(f"top_compose_pct = {top.frequency:.4f}")
    L(f"intensity = {r.compose_pairs[0].frequency if r.compose_pairs else 0:.4f}")
    L(f"compose_top5 = [")
    for cp in r.compose_pairs[:5]:
        L(f"  {cp.label}: {cp.frequency:.4f}")
    L(f"]")
    L(f"")

    pf = r.primary.frequencies
    mf = r.modifier.frequencies
    L(f"[LOGIC]")
    L(f"primary_dominant = {r.primary.dominant}")
    L(f"primary_pct = {r.primary.dominant_pct:.4f}")
    L(f"modifier_dominant = {r.modifier.dominant}")
    L(f"modifier_pct = {r.modifier.dominant_pct:.4f}")
    L(f"primary_FLOW = {pf.get('FLOW', 0):.4f}")
    L(f"primary_OBSERVE = {pf.get('OBSERVE', 0):.4f}")
    L(f"primary_LOGIC = {pf.get('LOGIC', 0):.4f}")
    L(f"primary_CHOOSE = {pf.get('CHOOSE', 0):.4f}")
    L(f"modifier_FLOW = {mf.get('FLOW', 0):.4f}")
    L(f"modifier_OBSERVE = {mf.get('OBSERVE', 0):.4f}")
    L(f"modifier_LOGIC = {mf.get('LOGIC', 0):.4f}")
    L(f"modifier_CHOOSE = {mf.get('CHOOSE', 0):.4f}")
    L(f"")

    L(f"[RUNTIME]")
    L(f"pattern = amino_acid_basis")
    L(f"aa_std = {r.aa_std:.4f}")
    L(f"aa_top5 = [")
    for aa, freq in list(r.aa_distribution.items())[:5]:
        props = AA_PROPERTIES.get(aa, {})
        L(f"  {aa}: {freq:.4f}  # {props.get('class', '?')}")
    L(f"]")
    L(f"")

    L(f"[MANIFEST]")
    L(f"format = {manifestation}")
    L(f"protein_length = {gene.protein_length}")
    L(f"stop_codon = {gene.stop_codon}")
    L(f"loss = {r.aa_std:.4f}")
    L(f"")

    sig_parts = [
        r.overall.dominant,
        r.wobble.dominant,
        obs_level,
        cycle_stab,
        dissolution,
        connection
    ]
    L(f"[SIGNATURE]")
    L(f"fingerprint = {'.'.join(sig_parts)}")

    return "\n".join(lines)


def write_dpl(result: TranslationResult, output_path: str):
    content = generate_dpl(result)
    with open(output_path, "w") as f:
        f.write(content)
