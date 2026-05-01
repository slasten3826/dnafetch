# fetch.py — FETCH layer
# Parses FASTA/GenBank, extracts CDS, validates sequences
# Part of DNAfetch v3

import gzip
import re
import os
from dataclasses import dataclass, field
from typing import List, Optional, Iterator, Tuple

# Standard genetic code
STOP_CODONS = {"TAA", "TAG", "TGA"}
START_CODON = "ATG"

CODON_TABLE = {
    "TTT": "Phe", "TTC": "Phe", "TTA": "Leu", "TTG": "Leu",
    "CTT": "Leu", "CTC": "Leu", "CTA": "Leu", "CTG": "Leu",
    "ATT": "Ile", "ATC": "Ile", "ATA": "Ile", "ATG": "Met",
    "GTT": "Val", "GTC": "Val", "GTA": "Val", "GTG": "Val",
    "TCT": "Ser", "TCC": "Ser", "TCA": "Ser", "TCG": "Ser",
    "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "TAT": "Tyr", "TAC": "Tyr", "TAA": "STOP", "TAG": "STOP",
    "CAT": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "AAT": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "GAT": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "TGT": "Cys", "TGC": "Cys", "TGA": "STOP", "TGG": "Trp",
    "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "AGT": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}

# Amino acid properties for functional layer
AA_PROPERTIES = {
    "Ala": {"class": "simple",    "hydrophobic": True,  "charge": "neutral"},
    "Arg": {"class": "basic",     "hydrophobic": False, "charge": "positive"},
    "Asn": {"class": "amide",     "hydrophobic": False, "charge": "neutral"},
    "Asp": {"class": "acidic",    "hydrophobic": False, "charge": "negative"},
    "Cys": {"class": "sulfur",    "hydrophobic": True,  "charge": "neutral"},
    "Gln": {"class": "amide",     "hydrophobic": False, "charge": "neutral"},
    "Glu": {"class": "acidic",    "hydrophobic": False, "charge": "negative"},
    "Gly": {"class": "simple",    "hydrophobic": True,  "charge": "neutral"},
    "His": {"class": "aromatic",  "hydrophobic": False, "charge": "positive"},
    "Ile": {"class": "aliphatic", "hydrophobic": True,  "charge": "neutral"},
    "Leu": {"class": "aliphatic", "hydrophobic": True,  "charge": "neutral"},
    "Lys": {"class": "basic",     "hydrophobic": False, "charge": "positive"},
    "Met": {"class": "sulfur",    "hydrophobic": True,  "charge": "neutral"},
    "Phe": {"class": "aromatic",  "hydrophobic": True,  "charge": "neutral"},
    "Pro": {"class": "cyclic",    "hydrophobic": True,  "charge": "neutral"},
    "Ser": {"class": "hydroxyl",  "hydrophobic": False, "charge": "neutral"},
    "Thr": {"class": "hydroxyl",  "hydrophobic": False, "charge": "neutral"},
    "Trp": {"class": "aromatic",  "hydrophobic": True,  "charge": "neutral"},
    "Tyr": {"class": "aromatic",  "hydrophobic": True,  "charge": "neutral"},
    "Val": {"class": "aliphatic", "hydrophobic": True,  "charge": "neutral"},
}


@dataclass
class Gene:
    """A single gene/transcript ready for ProcessLang translation."""
    name: str
    organism: str
    description: str
    sequence: str           # full mRNA/transcript
    cds: str                # coding sequence only
    codons: List[str]       # CDS split into triplets
    amino_acids: List[str]  # translated amino acids
    cds_start: int          # position of ATG in full sequence
    cds_end: int            # position of stop codon end
    transcript_id: str = ""
    gene_id: str = ""

    @property
    def cds_length(self) -> int:
        return len(self.cds)

    @property
    def codon_count(self) -> int:
        return len(self.codons)

    @property
    def protein_length(self) -> int:
        return len([aa for aa in self.amino_acids if aa != "STOP"])

    @property
    def stop_codon(self) -> str:
        if self.codons:
            last = self.codons[-1]
            if last in STOP_CODONS:
                return last
        return "none"


def parse_fasta_header(header: str) -> dict:
    """Extract gene name, organism, description from FASTA header."""
    info = {
        "name": "unknown",
        "organism": "unknown",
        "description": header.strip(">").strip(),
        "transcript_id": "",
        "gene_id": ""
    }

    # NM_XXXXXX.X format (RefSeq)
    refseq = re.match(r'>?(NM_\d+\.\d+)', header)
    if refseq:
        info["transcript_id"] = refseq.group(1)

    # Gene name in parentheses: (GRIN2A)
    gene_match = re.search(r'\(([A-Z0-9]+)\)', header)
    if gene_match:
        info["name"] = gene_match.group(1)

    # Organism: "Homo sapiens"
    org_match = re.search(r'(Homo sapiens|Mus musculus|Rattus norvegicus)', header)
    if org_match:
        info["organism"] = org_match.group(1)

    # GENCODE format: >ENST00000...| gene_name | ...
    enst_match = re.match(r'>?(ENST\d+\.\d+)\|', header)
    if enst_match:
        info["transcript_id"] = enst_match.group(1)
        parts = header.split("|")
        if len(parts) >= 6:
            info["gene_id"] = parts[1] if parts[1] else ""
            info["name"] = parts[5] if len(parts) > 5 and parts[5] else "unknown"

    return info


def find_longest_orf(sequence: str) -> Tuple[int, int]:
    """Find the longest open reading frame starting with ATG."""
    seq = sequence.upper().replace("U", "T")
    best_start = -1
    best_end = -1
    best_length = 0

    # Search all three reading frames
    for frame in range(3):
        i = frame
        while i < len(seq) - 2:
            codon = seq[i:i+3]
            if codon == START_CODON:
                # Found ATG, scan for stop
                j = i + 3
                while j < len(seq) - 2:
                    c = seq[j:j+3]
                    if c in STOP_CODONS:
                        orf_len = j + 3 - i
                        if orf_len > best_length:
                            best_length = orf_len
                            best_start = i
                            best_end = j + 3
                        break
                    j += 3
                # If no stop found, still check if longest
                else:
                    orf_len = j - i
                    if orf_len > best_length and orf_len > 300:  # min 100 aa
                        best_length = orf_len
                        best_start = i
                        best_end = j
            i += 3

    return best_start, best_end


def extract_cds(sequence: str) -> Tuple[str, int, int]:
    """Extract CDS from a transcript sequence."""
    seq = sequence.upper().replace("U", "T")

    start, end = find_longest_orf(seq)
    if start < 0:
        raise ValueError("No ORF found")

    cds = seq[start:end]
    return cds, start, end


def split_codons(cds: str) -> List[str]:
    """Split CDS into codons (triplets)."""
    return [cds[i:i+3] for i in range(0, len(cds), 3) if len(cds[i:i+3]) == 3]


def translate_codons(codons: List[str]) -> List[str]:
    """Translate codons to amino acids."""
    return [CODON_TABLE.get(c, "???") for c in codons]


def parse_gene(header: str, sequence: str) -> Gene:
    """Parse a single FASTA entry into a Gene object."""
    info = parse_fasta_header(header)
    cds, cds_start, cds_end = extract_cds(sequence)
    codons = split_codons(cds)
    amino_acids = translate_codons(codons)

    return Gene(
        name=info["name"],
        organism=info["organism"],
        description=info["description"],
        sequence=sequence.upper().replace("U", "T"),
        cds=cds,
        codons=codons,
        amino_acids=amino_acids,
        cds_start=cds_start,
        cds_end=cds_end,
        transcript_id=info["transcript_id"],
        gene_id=info["gene_id"]
    )


def read_fasta(filepath: str) -> Iterator[Gene]:
    """Read FASTA file (plain or gzipped), yield Gene objects."""
    opener = gzip.open if filepath.endswith(".gz") else open
    mode = "rt" if filepath.endswith(".gz") else "r"

    header = None
    seq_lines = []

    with opener(filepath, mode) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    sequence = "".join(seq_lines)
                    try:
                        yield parse_gene(header, sequence)
                    except ValueError:
                        pass  # skip entries with no ORF
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)

        # last entry
        if header is not None:
            sequence = "".join(seq_lines)
            try:
                yield parse_gene(header, sequence)
            except ValueError:
                pass


def fetch_single(filepath: str) -> Gene:
    """Fetch a single gene from a FASTA file."""
    for gene in read_fasta(filepath):
        return gene
    raise ValueError(f"No valid gene found in {filepath}")


def fetch_batch(filepath: str, limit: int = 0) -> Iterator[Gene]:
    """Fetch all genes from a FASTA file (supports .gz)."""
    count = 0
    for gene in read_fasta(filepath):
        yield gene
        count += 1
        if limit > 0 and count >= limit:
            break
