const THREE_TO_ONE = {
  ALA: "A", ARG: "R", ASN: "N", ASP: "D", CYS: "C", GLN: "Q", GLU: "E",
  GLY: "G", HIS: "H", ILE: "I", LEU: "L", LYS: "K", MET: "M", PHE: "F",
  PRO: "P", SER: "S", THR: "T", TRP: "W", TYR: "Y", VAL: "V",
  MSE: "M", SEC: "U", PYL: "O",
};

const THREE_TO_FULL_NAME = {
  ALA: "Alanine", ARG: "Arginine", ASN: "Asparagine", ASP: "Aspartate",
  CYS: "Cysteine", GLN: "Glutamine", GLU: "Glutamate", GLY: "Glycine",
  HIS: "Histidine", ILE: "Isoleucine", LEU: "Leucine", LYS: "Lysine",
  MET: "Methionine", PHE: "Phenylalanine", PRO: "Proline", SER: "Serine",
  THR: "Threonine", TRP: "Tryptophan", TYR: "Tyrosine", VAL: "Valine",
  MSE: "Selenomethionine", SEC: "Selenocysteine", PYL: "Pyrrolysine",
};

let viewer = null;
let proteins = [];
let properties = {};
let liveMeta = {};
let currentStyle = "stick";
let currentPdbId = null;
let residues = [];
let highlightedResidue = null;
let searchDebounceTimer = null;
let searchRequestSeq = 0;

async function init() {
  viewer = $3Dmol.createViewer(document.getElementById("viewer"), {
    backgroundColor: "#0f1420",
  });

  const [proteinsRes, propertiesRes] = await Promise.all([
    fetch("/api/proteins"),
    fetch("/api/properties"),
  ]);
  proteins = await proteinsRes.json();
  properties = await propertiesRes.json();

  renderProteinList();
  document.getElementById("style-select").addEventListener("change", (e) => {
    currentStyle = e.target.value;
    applyStyle();
  });
  document.getElementById("search-input").addEventListener("input", onSearchInput);
  document.getElementById("rotate-btn").addEventListener("click", rotate360);

  if (proteins.length > 0) {
    selectProtein(proteins[0].pdb_id);
  }
}

function rotate360() {
  if (!viewer) return;
  const btn = document.getElementById("rotate-btn");
  const durationMs = 3600;
  const totalDegrees = 360;
  let lastTime = null;
  let rotated = 0;

  btn.disabled = true;

  // drive rotation manually via requestAnimationFrame so every step is
  // guaranteed to repaint - 3Dmol's own rotate()/spin() timers can silently
  // skip intermediate repaints
  function step(now) {
    if (lastTime === null) lastTime = now;
    const delta = (totalDegrees * (now - lastTime)) / durationMs;
    lastTime = now;
    rotated += delta;

    viewer.rotate(delta, "y");
    viewer.render();

    if (rotated < totalDegrees) {
      requestAnimationFrame(step);
    } else {
      btn.disabled = false;
    }
  }
  requestAnimationFrame(step);
}

function onSearchInput(e) {
  const query = e.target.value.trim();
  clearTimeout(searchDebounceTimer);

  if (query.length < 2) {
    document.getElementById("search-results").innerHTML = "";
    document.getElementById("search-status").textContent = "";
    return;
  }

  document.getElementById("search-status").textContent = "Searching\u2026";
  searchDebounceTimer = setTimeout(() => runSearch(query), 350);
}

async function runSearch(query) {
  const requestId = ++searchRequestSeq;
  const status = document.getElementById("search-status");
  try {
    const response = await fetch(`/api/search?q=${encodeURIComponent(query)}`);
    const results = await response.json();
    if (requestId !== searchRequestSeq) return; // a newer search superseded this one

    status.textContent = results.length
      ? `${results.length} result${results.length === 1 ? "" : "s"}`
      : "No matches";
    renderSearchResults(results);
  } catch (err) {
    if (requestId !== searchRequestSeq) return;
    status.textContent = "Search failed";
  }
}

function renderSearchResults(results) {
  const list = document.getElementById("search-results");
  list.innerHTML = "";
  results.forEach((r) => {
    const li = document.createElement("li");
    li.dataset.pdbId = r.pdb_id;

    const title = document.createElement("span");
    title.className = "result-title";
    title.textContent = r.title || r.pdb_id;
    li.appendChild(title);

    const idLine = document.createElement("span");
    idLine.className = "result-id";
    idLine.textContent = r.pdb_id;
    li.appendChild(idLine);

    li.addEventListener("click", () => selectLiveProtein(r.pdb_id, r));
    list.appendChild(li);
  });
}

function renderProteinList() {
  const list = document.getElementById("protein-list");
  list.innerHTML = "";
  proteins.forEach((p) => {
    const li = document.createElement("li");
    li.textContent = `${p.name} (${p.pdb_id})`;
    li.dataset.pdbId = p.pdb_id;
    li.addEventListener("click", () => selectProtein(p.pdb_id));
    list.appendChild(li);
  });
}

function markActive(pdbId) {
  document.querySelectorAll("#protein-list li").forEach((li) => {
    li.classList.toggle("active", li.dataset.pdbId === pdbId);
  });
  document.querySelectorAll("#search-results li").forEach((li) => {
    li.classList.toggle("active", li.dataset.pdbId === pdbId);
  });
}

async function selectProtein(pdbId) {
  const response = await fetch(`/data/structures/${pdbId}.pdb`);
  const pdbText = await response.text();
  loadModel(pdbId, pdbText);
  renderProperties(pdbId);
}

async function selectLiveProtein(pdbId, resultMeta) {
  const status = document.getElementById("search-status");
  status.textContent = `Loading ${pdbId}\u2026`;
  try {
    const [structureRes, infoRes] = await Promise.all([
      fetch(`/api/live/${pdbId}/structure`),
      fetch(`/api/live/${pdbId}/info`),
    ]);
    if (!structureRes.ok || !infoRes.ok) throw new Error("fetch failed");

    const format = structureRes.headers.get("X-Structure-Format") || "pdb";
    const structureText = await structureRes.text();
    const info = await infoRes.json();

    liveMeta[pdbId] = {
      pdb_id: pdbId,
      name: resultMeta.title || pdbId,
      title: info.title,
      experimental_method: info.experimental_method,
      resolution_angstrom: info.resolution_angstrom,
    };
    properties[pdbId] = info.properties;

    loadModel(pdbId, structureText, format);
    renderProperties(pdbId);
    status.textContent = `Loaded ${pdbId}`;
  } catch (err) {
    status.textContent = `Could not load ${pdbId}`;
  }
}

function loadModel(pdbId, structureText, format = "pdb") {
  currentPdbId = pdbId;
  markActive(pdbId);
  document.getElementById("viewer-pdbid-badge").textContent = pdbId;

  viewer.clear();
  const model = viewer.addModel(structureText, format);
  highlightedResidue = null;
  applyStyle();
  viewer.zoomTo();
  viewer.render();

  buildSequencePanel(model);
}

function buildSequencePanel(model) {
  const seen = new Set();
  residues = [];
  model.selectedAtoms({}).forEach((atom) => {
    if (atom.hetflag) return; // skip waters/ligands, keep only chain residues
    const key = `${atom.chain}_${atom.resi}`;
    if (seen.has(key)) return;
    seen.add(key);
    residues.push({ chain: atom.chain, resi: atom.resi, resn: atom.resn });
  });
  residues.sort((a, b) =>
    a.chain !== b.chain ? a.chain.localeCompare(b.chain) : a.resi - b.resi
  );
  renderSequenceStrip();
}

function renderSequenceStrip() {
  const container = document.getElementById("sequence-strip");
  container.innerHTML = "";
  let lastChain = null;
  residues.forEach((res) => {
    if (res.chain !== lastChain) {
      lastChain = res.chain;
      const label = document.createElement("div");
      label.className = "chain-label";
      label.textContent = `Chain ${res.chain}`;
      container.appendChild(label);
    }
    const span = document.createElement("span");
    span.className = "residue";
    span.textContent = THREE_TO_ONE[res.resn] || "X";
    span.title = `${res.resn} ${res.resi} (chain ${res.chain})`;
    span.addEventListener("click", () => selectResidue(res.chain, res.resi, res.resn, span));
    container.appendChild(span);
  });
}

function selectResidue(chain, resi, resn, spanEl) {
  document
    .querySelectorAll("#sequence-strip .residue.active")
    .forEach((el) => el.classList.remove("active"));
  spanEl.classList.add("active");
  spanEl.scrollIntoView({ block: "nearest", inline: "nearest" });

  highlightedResidue = { chain, resi };
  applyHighlight();
  viewer.render();
  viewer.center({ chain, resi }, 500);

  showResidueTooltip(spanEl, resn, resi, chain);
}

function showResidueTooltip(spanEl, resn, resi, chain) {
  const tooltip = document.getElementById("residue-tooltip");
  const name = THREE_TO_FULL_NAME[resn] || resn;
  tooltip.textContent = `${name} (${resn}) \u2014 residue ${resi}, chain ${chain}`;

  const rect = spanEl.getBoundingClientRect();
  tooltip.style.left = `${rect.left + rect.width / 2}px`;
  tooltip.style.top = `${rect.top - 8}px`;
  tooltip.hidden = false;

  clearTimeout(showResidueTooltip.timer);
  showResidueTooltip.timer = setTimeout(() => {
    tooltip.hidden = true;
  }, 2500);
}

function applyHighlight() {
  if (!viewer || !highlightedResidue) return;
  const { chain, resi } = highlightedResidue;
  viewer.setStyle(
    { chain, resi },
    {
      stick: { radius: 0.35, colorscheme: "Jmol" },
      sphere: { scale: 0.45, colorscheme: "Jmol" },
    }
  );
}

function applyStyle() {
  if (!viewer) return;
  viewer.setStyle({}, {});
  const styleMap = {
    stick: { stick: { radius: 0.15 } },
    sphere: { sphere: { scale: 0.3 } },
    cartoon: { cartoon: { color: "spectrum" } },
    line: { line: {} },
  };
  // hide crystallographic waters (HOH), which are lone unbonded oxygens
  viewer.setStyle(
    { resn: "HOH", invert: true },
    styleMap[currentStyle] || styleMap.stick
  );
  applyHighlight();
  viewer.render();
}

function renderProperties(pdbId) {
  const meta = proteins.find((p) => p.pdb_id === pdbId) || liveMeta[pdbId];
  const props = properties[pdbId];

  document.getElementById("properties-title").textContent = meta
    ? `${meta.name} (${pdbId})`
    : pdbId;
  document.getElementById("properties-description").textContent = meta
    ? meta.title || meta.description
    : "";

  const dl = document.getElementById("properties-list");
  dl.innerHTML = "";

  if (!props) {
    dl.innerHTML = "<dd>No computed properties available.</dd>";
    return;
  }

  const rows = [
    ["Sequence length", `${props.sequence_length} aa`],
    ["Molecular weight", `${props.molecular_weight_da.toLocaleString()} Da`],
    ["Isoelectric point (pI)", props.isoelectric_point],
    ["Net charge @ pH 7", props.net_charge_at_ph7],
    ["GRAVY (hydropathy)", props.gravy_hydropathy],
    ["Aromaticity", props.aromaticity],
    ["Instability index", props.instability_index],
    [
      "Stability prediction",
      props.is_stable_prediction ? "Stable" : "Potentially unstable",
    ],
    [
      "Secondary structure",
      `Helix ${(props.secondary_structure_fraction.helix * 100).toFixed(1)}% / ` +
        `Sheet ${(props.secondary_structure_fraction.sheet * 100).toFixed(1)}% / ` +
        `Turn ${(props.secondary_structure_fraction.turn * 100).toFixed(1)}%`,
    ],
    [
      "Extinction coeff. (reduced / oxidized)",
      `${props.molar_extinction_coefficient.reduced_cysteines} / ` +
        `${props.molar_extinction_coefficient.oxidized_cystines} M\u207B\u00B9cm\u207B\u00B9`,
    ],
  ];

  if (meta && meta.experimental_method) {
    rows.push(["Experimental method", meta.experimental_method]);
  }
  if (meta && meta.resolution_angstrom) {
    rows.push(["Resolution", `${meta.resolution_angstrom} \u00C5`]);
  }

  rows.forEach(([label, value]) => {
    const dt = document.createElement("dt");
    dt.textContent = label;
    const dd = document.createElement("dd");
    dd.textContent = value;
    dl.appendChild(dt);
    dl.appendChild(dd);
  });
}

init();
