"""Compute predicted physicochemical properties for each protein's sequence
using Biopython's ProtParam module.

Usage:
    python scripts/compute_properties.py
"""

import json
import re
from pathlib import Path

from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

DATA_DIR = Path(__file__).parent.parent / "data"
SEQUENCES_DIR = DATA_DIR / "sequences"

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def longest_sequence(fasta_path: Path) -> str:
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        return ""
    longest = max(records, key=lambda r: len(r.seq))
    return str(longest.seq)


def clean_sequence(sequence: str) -> str:
    """Strip residues ProtParam can't handle (non-standard/ambiguous codes)."""
    return re.sub(f"[^{''.join(VALID_AA)}]", "", sequence.upper())


def analyze_sequence(sequence: str) -> dict:
    analysis = ProteinAnalysis(sequence)
    helix, turn, sheet = analysis.secondary_structure_fraction()
    instability = analysis.instability_index()
    reduced, oxidized = analysis.molar_extinction_coefficient()

    return {
        "sequence_length": len(sequence),
        "molecular_weight_da": round(analysis.molecular_weight(), 1),
        "isoelectric_point": round(analysis.isoelectric_point(), 2),
        "aromaticity": round(analysis.aromaticity(), 4),
        "instability_index": round(instability, 2),
        "is_stable_prediction": instability < 40,
        "gravy_hydropathy": round(analysis.gravy(), 4),
        "net_charge_at_ph7": round(analysis.charge_at_pH(7.0), 2),
        "secondary_structure_fraction": {
            "helix": round(helix, 3),
            "turn": round(turn, 3),
            "sheet": round(sheet, 3),
        },
        "molar_extinction_coefficient": {
            "reduced_cysteines": reduced,
            "oxidized_cystines": oxidized,
        },
        "amino_acid_percent": {
            aa: round(pct / 100, 4)
            for aa, pct in analysis.amino_acids_percent.items()
        },
    }


def main() -> None:
    proteins = json.loads((DATA_DIR / "proteins.json").read_text())
    properties = {}

    for entry in proteins:
        pdb_id = entry["pdb_id"]
        fasta_path = SEQUENCES_DIR / f"{pdb_id}.fasta"
        if not fasta_path.exists():
            print(f"! no sequence file for {pdb_id}, skipping")
            continue

        raw_sequence = longest_sequence(fasta_path)
        sequence = clean_sequence(raw_sequence)
        if not sequence:
            print(f"! no analyzable sequence for {pdb_id}, skipping")
            continue

        print(f"Analyzing {pdb_id} ({entry['name']}), {len(sequence)} residues...")
        properties[pdb_id] = analyze_sequence(sequence)

    output_path = DATA_DIR / "properties.json"
    output_path.write_text(json.dumps(properties, indent=2))
    print(f"\nSaved computed properties to {output_path}")


if __name__ == "__main__":
    main()
