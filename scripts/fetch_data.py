"""Download structure files, sequences, and metadata for the project's 10 proteins
from the RCSB Protein Data Bank (open data source, no API key required).

Usage:
    python scripts/fetch_data.py
"""

import json
import sys
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).parent))
from proteins import PROTEINS

DATA_DIR = Path(__file__).parent.parent / "data"
STRUCTURES_DIR = DATA_DIR / "structures"
SEQUENCES_DIR = DATA_DIR / "sequences"
METADATA_DIR = DATA_DIR / "metadata"

STRUCTURE_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
FASTA_URL = "https://www.rcsb.org/fasta/entry/{pdb_id}"
ENTRY_METADATA_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"


def download(url: str, dest: Path, timeout: int = 30) -> bool:
    try:
        response = requests.get(url, timeout=timeout)
        response.raise_for_status()
    except requests.RequestException as exc:
        print(f"  ! failed to download {url}: {exc}")
        return False
    dest.write_bytes(response.content)
    return True


def fetch_protein(entry: dict) -> dict:
    pdb_id = entry["pdb_id"]
    print(f"Fetching {pdb_id} ({entry['name']})...")

    structure_path = STRUCTURES_DIR / f"{pdb_id}.pdb"
    download(STRUCTURE_URL.format(pdb_id=pdb_id), structure_path)

    sequence_path = SEQUENCES_DIR / f"{pdb_id}.fasta"
    download(FASTA_URL.format(pdb_id=pdb_id), sequence_path)

    title = None
    resolution = None
    method = None
    try:
        response = requests.get(ENTRY_METADATA_URL.format(pdb_id=pdb_id), timeout=30)
        response.raise_for_status()
        meta = response.json()
        (METADATA_DIR / f"{pdb_id}.json").write_text(json.dumps(meta, indent=2))
        title = meta.get("struct", {}).get("title")
        resolution_list = meta.get("rcsb_entry_info", {}).get("resolution_combined")
        resolution = resolution_list[0] if resolution_list else None
        method_list = meta.get("exptl", [])
        method = method_list[0].get("method") if method_list else None
    except requests.RequestException as exc:
        print(f"  ! failed to fetch metadata for {pdb_id}: {exc}")

    return {
        **entry,
        "title": title,
        "resolution_angstrom": resolution,
        "experimental_method": method,
    }


def main() -> None:
    for directory in (STRUCTURES_DIR, SEQUENCES_DIR, METADATA_DIR):
        directory.mkdir(parents=True, exist_ok=True)

    enriched = [fetch_protein(entry) for entry in PROTEINS]

    output_path = DATA_DIR / "proteins.json"
    output_path.write_text(json.dumps(enriched, indent=2))
    print(f"\nSaved combined protein index to {output_path}")


if __name__ == "__main__":
    main()
