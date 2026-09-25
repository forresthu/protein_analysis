"""Flask app serving the protein 3D viewer and property dashboard."""

import io
import json
import re
from pathlib import Path

import requests
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from flask import Flask, Response, abort, jsonify, render_template, request, send_from_directory

BASE_DIR = Path(__file__).parent
DATA_DIR = BASE_DIR / "data"
STRUCTURES_DIR = DATA_DIR / "structures"
SEQUENCES_DIR = DATA_DIR / "sequences"
for _dir in (STRUCTURES_DIR, SEQUENCES_DIR):
    _dir.mkdir(parents=True, exist_ok=True)

PDB_ID_RE = re.compile(r"^[A-Za-z0-9]{4}$")
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

STRUCTURE_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
CIF_URL = "https://files.rcsb.org/download/{pdb_id}.cif"
FASTA_URL = "https://www.rcsb.org/fasta/entry/{pdb_id}"
ENTRY_METADATA_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL_URL = "https://data.rcsb.org/graphql"

app = Flask(__name__)


def load_json(name: str):
    path = DATA_DIR / name
    if not path.exists():
        return None
    return json.loads(path.read_text())


def validate_pdb_id(pdb_id: str) -> str:
    if not PDB_ID_RE.match(pdb_id or ""):
        abort(400, "invalid PDB id")
    return pdb_id.upper()


def fetch_structure_text(pdb_id: str) -> tuple[str, str]:
    """Returns (structure_text, format) where format is 'pdb' or 'cif'.

    Some large structures (e.g. cryo-EM assemblies) exceed the legacy PDB
    format's limits and are only distributed as mmCIF, so we fall back to it.
    """
    pdb_path = STRUCTURES_DIR / f"{pdb_id}.pdb"
    cif_path = STRUCTURES_DIR / f"{pdb_id}.cif"
    if pdb_path.exists():
        return pdb_path.read_text(), "pdb"
    if cif_path.exists():
        return cif_path.read_text(), "cif"

    response = requests.get(STRUCTURE_URL.format(pdb_id=pdb_id), timeout=15)
    if response.status_code == 200:
        pdb_path.write_text(response.text)
        return response.text, "pdb"

    response = requests.get(CIF_URL.format(pdb_id=pdb_id), timeout=15)
    response.raise_for_status()
    cif_path.write_text(response.text)
    return response.text, "cif"


def fetch_sequence_text(pdb_id: str) -> str:
    path = SEQUENCES_DIR / f"{pdb_id}.fasta"
    if path.exists():
        return path.read_text()
    response = requests.get(FASTA_URL.format(pdb_id=pdb_id), timeout=15)
    response.raise_for_status()
    path.write_text(response.text)
    return response.text


def clean_sequence(sequence: str) -> str:
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
            aa: round(pct / 100, 4) for aa, pct in analysis.amino_acids_percent.items()
        },
    }


def compute_properties_for(pdb_id: str):
    cache = load_json("properties.json") or {}
    if pdb_id in cache:
        return cache[pdb_id]

    fasta_text = fetch_sequence_text(pdb_id)
    records = list(SeqIO.parse(io.StringIO(fasta_text), "fasta"))
    if not records:
        return None
    longest = max(records, key=lambda r: len(r.seq))
    sequence = clean_sequence(str(longest.seq))
    if not sequence:
        return None

    props = analyze_sequence(sequence)
    cache[pdb_id] = props
    (DATA_DIR / "properties.json").write_text(json.dumps(cache, indent=2))
    return props


def fetch_entry_meta(pdb_id: str) -> dict:
    response = requests.get(ENTRY_METADATA_URL.format(pdb_id=pdb_id), timeout=15)
    response.raise_for_status()
    meta = response.json()
    resolution_list = meta.get("rcsb_entry_info", {}).get("resolution_combined")
    method_list = meta.get("exptl", [])
    return {
        "title": meta.get("struct", {}).get("title"),
        "resolution_angstrom": resolution_list[0] if resolution_list else None,
        "experimental_method": method_list[0].get("method") if method_list else None,
    }


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/api/proteins")
def api_proteins():
    proteins = load_json("proteins.json") or []
    return jsonify(proteins)


@app.route("/api/properties")
def api_properties():
    properties = load_json("properties.json") or {}
    return jsonify(properties)


@app.route("/api/search")
def api_search():
    query = request.args.get("q", "").strip()[:100]
    if not query:
        return jsonify([])

    search_payload = {
        "query": {"type": "terminal", "service": "full_text", "parameters": {"value": query}},
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": 8}},
    }
    try:
        search_resp = requests.post(SEARCH_URL, json=search_payload, timeout=10)
        search_resp.raise_for_status()
        ids = [hit["identifier"] for hit in search_resp.json().get("result_set", [])]
    except requests.RequestException:
        return jsonify([])
    if not ids:
        return jsonify([])

    graphql_query = (
        "{ entries(entry_ids: %s) { rcsb_id struct { title } exptl { method } "
        "rcsb_entry_info { resolution_combined } } }" % json.dumps(ids)
    )
    try:
        meta_resp = requests.post(GRAPHQL_URL, json={"query": graphql_query}, timeout=10)
        meta_resp.raise_for_status()
        entries = meta_resp.json().get("data", {}).get("entries") or []
    except requests.RequestException:
        entries = []

    results = []
    for entry in entries:
        resolution_list = (entry.get("rcsb_entry_info") or {}).get("resolution_combined")
        method_list = entry.get("exptl") or []
        results.append(
            {
                "pdb_id": entry["rcsb_id"],
                "title": (entry.get("struct") or {}).get("title"),
                "resolution_angstrom": resolution_list[0] if resolution_list else None,
                "experimental_method": method_list[0].get("method") if method_list else None,
            }
        )
    order = {pid: i for i, pid in enumerate(ids)}
    results.sort(key=lambda r: order.get(r["pdb_id"], len(ids)))
    return jsonify(results)


@app.route("/api/live/<pdb_id>/structure")
def api_live_structure(pdb_id):
    pdb_id = validate_pdb_id(pdb_id)
    try:
        text, fmt = fetch_structure_text(pdb_id)
    except requests.RequestException:
        abort(404)
    response = Response(text, mimetype="text/plain")
    response.headers["X-Structure-Format"] = fmt
    return response


@app.route("/api/live/<pdb_id>/info")
def api_live_info(pdb_id):
    pdb_id = validate_pdb_id(pdb_id)
    try:
        meta = fetch_entry_meta(pdb_id)
        properties = compute_properties_for(pdb_id)
    except requests.RequestException:
        abort(404)
    if properties is None:
        abort(404)
    return jsonify({**meta, "pdb_id": pdb_id, "properties": properties})


@app.route("/data/structures/<path:filename>")
def structure_file(filename):
    return send_from_directory(DATA_DIR / "structures", filename)


if __name__ == "__main__":
    # avoid 5000: macOS AirPlay Receiver/ControlCenter also listens there and
    # can intercept requests, returning 403 Forbidden instead of reaching Flask.
    app.run(debug=True, port=8000)
