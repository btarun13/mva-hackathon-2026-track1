"""Tier B: HGNC mapping, protein-changing variants, STRING reference-vs-proband interactome, BUB1B domain map,
DGIdb + ChEMBL druggability. Reads work/vep.tsv (panel VEP table); writes results/.
Run: uv run --no-project --with pandas,requests,networkx,matplotlib python scripts/tier_b_interactome.py
"""
import json, sys, time
from pathlib import Path
import pandas as pd, requests, networkx as nx
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from mva_track1 import read_vep_tab
R = ROOT / "results"; R.mkdir(exist_ok=True)
S = requests.Session(); S.headers["User-Agent"] = "mva-hackathon-tier-b"

def get(url, **kw):
    for i in range(3):
        try:
            r = S.get(url, timeout=60, **kw); r.raise_for_status(); return r
        except Exception as e:
            if i == 2: raise
            time.sleep(2)

# ---------------------------------------------------------------- 1. HGNC mapping of every gene in vep.tsv
vep = read_vep_tab(str(ROOT / "work/vep.tsv"))
gt = pd.read_csv(ROOT / "work/subset.gt.tsv", sep="\t", header=None,
                 names=["ID","CHROM","POS","REF","ALT","FILTER","GT","DP","GQ","AD"], dtype=str)
vep = vep.merge(gt, left_on="Uploaded_variation", right_on="ID", how="left")
rows = []
for ensg in sorted(vep["Gene"].dropna().unique()):
    d = get(f"https://rest.genenames.org/fetch/ensembl_gene_id/{ensg}", headers={"Accept": "application/json"}).json()
    docs = d["response"]["docs"]
    if docs:
        x = docs[0]
        rows.append(dict(ensembl_gene_id=ensg, hgnc_id=x.get("hgnc_id"), symbol=x.get("symbol"), name=x.get("name"),
                         uniprot=";".join(x.get("uniprot_ids", []) or []), locus_type=x.get("locus_type"),
                         alias=";".join(x.get("alias_symbol", []) or []), location=x.get("location")))
    else:
        rows.append(dict(ensembl_gene_id=ensg, hgnc_id="", symbol="", name="(not in HGNC)", uniprot="", locus_type="", alias="", location=""))
hgnc = pd.DataFrame(rows); hgnc.to_csv(R / "hgnc_map.tsv", sep="\t", index=False)
print("HGNC map:", len(hgnc), "genes;", (hgnc.symbol == "").sum(), "without a symbol")

# ---------------------------------------------------------------- 2. Protein-changing variants in the proband
PROT = {"missense_variant","stop_gained","stop_lost","start_lost","frameshift_variant","inframe_insertion",
        "inframe_deletion","protein_altering_variant","splice_acceptor_variant","splice_donor_variant"}
vep["max_af"] = pd.to_numeric(vep["MAX_AF"], errors="coerce")
pc = vep[vep["Consequence"].apply(lambda c: bool(set(str(c).split(",")) & PROT))].copy()
pc["rare_lt_1pct"] = pc["max_af"].fillna(0) < 0.01
pc = pc[["Uploaded_variation","SYMBOL","Gene","Consequence","IMPACT","Protein_position","Amino_acids","HGVSc","HGVSp",
         "am_class","am_pathogenicity","ClinVar_CLNSIG","max_af","GT","rare_lt_1pct"]].sort_values("max_af", na_position="first")
pc.to_csv(R / "protein_changing_variants.tsv", sep="\t", index=False)
altered = sorted(pc.loc[pc.rare_lt_1pct, "SYMBOL"].dropna().unique())
print("protein-changing variants:", len(pc), "| rare (<1%):", int(pc.rare_lt_1pct.sum()), "| altered proteins:", altered)

# ---------------------------------------------------------------- 3. STRING network: reference vs proband
SEED = ["BUB1B","BUB1","BUB3","MAD1L1","MAD2L1","CDC20","TRIP13","MAD2L1BP","CEP57","CENPE","KIF18A","CENATAC",
        "TTK","PLK1","AURKB","KNL1","PPP2R5B","PPP2R5C","PPP2R5E","PPP2CA","CDC27","ANAPC1","CDK1","ZW10","CASC5"]
SEED = [g for g in SEED if g != "CASC5"]  # CASC5 = KNL1 old symbol
SP, SCORE = 9606, 700
ids = get("https://string-db.org/api/tsv/get_string_ids", params={"identifiers": "\r".join(SEED), "species": SP, "limit": 1}).text
ids = pd.read_csv(pd.io.common.StringIO(ids), sep="\t")
sid = dict(zip(ids.queryIndex.map(lambda i: SEED[i]), ids.stringId))
net = get("https://string-db.org/api/tsv/network", params={"identifiers": "\r".join(sid.values()), "species": SP,
          "required_score": SCORE, "network_type": "functional"}).text
net = pd.read_csv(pd.io.common.StringIO(net), sep="\t").drop_duplicates(subset=["preferredName_A","preferredName_B"])
net = net[["preferredName_A","preferredName_B","score","escore","dscore","tscore","ascore","pscore","nscore","fscore"]]
net.to_csv(R / "string_edges_reference.tsv", sep="\t", index=False)
enr = get("https://string-db.org/api/tsv/enrichment", params={"identifiers": "\r".join(sid.values()), "species": SP}).text
enr = pd.read_csv(pd.io.common.StringIO(enr), sep="\t")
enr[enr.fdr < 0.05].sort_values("fdr").to_csv(R / "string_enrichment.tsv", sep="\t", index=False)
print("STRING: nodes", len(sid), "edges", len(net), "| enriched terms (FDR<0.05):", int((enr.fdr < 0.05).sum()))

# ---------------------------------------------------------------- 4. Interactome diff for altered proteins (BUB1B)
# BUBR1 (O60566, 1050 aa) regions and the partners each one mediates (UniProt + literature: Overlack 2015, Suijkerbuijk 2012,
# Kruse 2013, Di Fiore 2015, Chao 2012, Lara-Gonzalez 2011, Mao 2003/ Ciossani 2018 for CENPE).
BUB1B_REGIONS = [
    ("KEN1 box",            26,   28, ["CDC20","MAD2L1"], "MCC assembly; CDC20 pseudosubstrate"),
    ("TPR domain",          62,  226, ["KNL1","BUB1"],    "kinetochore recruitment via KNL1 MELT/KI motifs"),
    ("D-box",              224,  232, ["CDC20","CDC27","ANAPC1"], "APC/C-CDC20 inhibition"),
    ("KEN2 box",           304,  306, ["CDC20"],          "second CDC20 site; APC/C inhibition"),
    ("GLEBS motif",        400,  431, ["BUB3"],           "BUB3 binding; kinetochore localisation"),
    ("ABBA motifs",        525,  577, ["CDC20"],          "CDC20 binding (ABBA1/2)"),
    ("KARD",               665,  682, ["PPP2R5B","PPP2R5C","PPP2R5E","PPP2CA","PLK1","CDK1","TTK"], "PP2A-B56 recruitment; phospho-sites for CDK1/PLK1/MPS1"),
    ("Protein kinase (pseudokinase)", 766, 1050, ["CENPE"], "CENPE binding; required for BUBR1 stability"),
]
TRUNC_END = 736   # p.Leu737Ter: residues 1-736 retained if the transcript escapes NMD
MIS_POS = 1002    # p.Asn1002Lys
part2reg = {}
for name, a, b, partners, fn in BUB1B_REGIONS:
    for p in partners: part2reg.setdefault(p, []).append((name, a, b, fn))

def classify(partner):
    regs = part2reg.get(partner)
    if not regs:
        return "weakened (dosage)", "region not annotated", "allele 1: no protein (NMD) | allele 2: destabilised protein"
    names = ";".join(r[0] for r in regs); fn = ";".join(r[3] for r in regs)
    lost_on_trunc = any(r[1] > TRUNC_END for r in regs)
    hit_by_mis = any(r[1] <= MIS_POS <= r[2] for r in regs)
    if lost_on_trunc or hit_by_mis:
        return "lost / severely impaired", names, fn
    return "weakened (dosage; binding region retained)", names, fn

diff = []
for _, e in net.iterrows():
    a, b = e.preferredName_A, e.preferredName_B
    for alt in altered:
        if alt in (a, b):
            partner = b if a == alt else a
            state, region, fn = classify(partner)
            diff.append(dict(altered_protein=alt, partner=partner, string_score=e.score, experimental=e.escore, database=e.dscore,
                             textmining=e.tscore, mediating_region=region, function=fn, predicted_state_proband=state))
diff = pd.DataFrame(diff).sort_values(["predicted_state_proband","string_score"], ascending=[True, False])
diff.to_csv(R / "interactome_diff.tsv", sep="\t", index=False)
print("interactome diff rows:", len(diff)); print(diff.predicted_state_proband.value_counts().to_string())

# ---------------------------------------------------------------- 5. Figures: reference vs proband network, domain map
G = nx.Graph()
for _, e in net.iterrows(): G.add_edge(e.preferredName_A, e.preferredName_B, w=e.score)
pos = nx.spring_layout(G, seed=7, k=0.9)
state = {(r.altered_protein, r.partner): r.predicted_state_proband for r in diff.itertuples()}
fig, axes = plt.subplots(1, 2, figsize=(18, 9))
for ax, title, proband in zip(axes, ["Reference interactome (no protein-changing variant)", "Proband (BUB1B p.Leu737Ter / p.Asn1002Lys)"], [False, True]):
    for u, v, d in G.edges(data=True):
        st = state.get((u, v)) or state.get((v, u))
        if proband and st:
            style = ":" if st.startswith("lost") else "--"; col = "#c0392b" if st.startswith("lost") else "#e67e22"; lw = 1.2
        else:
            style, col, lw = "-", "#95a5a6", 0.5 + 2.5 * (d["w"] - 0.7) / 0.3
        ax.plot(*zip(pos[u], pos[v]), linestyle=style, color=col, linewidth=lw, zorder=1)
    cols = ["#c0392b" if (proband and n in altered) else ("#2980b9" if n in {"BUB1B","CEP57","TRIP13"} else "#bdc3c7") for n in G.nodes]
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=cols, node_size=900, edgecolors="k", linewidths=0.6)
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=8)
    ax.set_title(title); ax.axis("off")
axes[1].plot([], [], ":", color="#c0392b", label="lost / severely impaired"); axes[1].plot([], [], "--", color="#e67e22", label="weakened (dosage)")
axes[1].legend(loc="lower left", fontsize=8)
fig.suptitle(f"STRING v12 functional network, score ≥ {SCORE/1000:.1f}; {len(G)} proteins, {G.number_of_edges()} edges", y=0.98)
fig.tight_layout(); fig.savefig(R / "interactome_reference_vs_proband.png", dpi=150); plt.close(fig)

fig, ax = plt.subplots(figsize=(14, 3.2))
ax.add_patch(plt.Rectangle((1, 0.4), 1050, 0.2, color="#dfe6e9"))
palette = ["#74b9ff","#55efc4","#ffeaa7","#fab1a0","#a29bfe","#fd79a8","#81ecec","#636e72"]
for (name, a, b, partners, fn), c in zip(BUB1B_REGIONS, palette):
    ax.add_patch(plt.Rectangle((a, 0.35), max(b - a, 6), 0.3, color=c, ec="k", lw=0.4))
    ax.text((a + b) / 2, 0.72 if len(name) < 10 else 0.80, name, ha="center", fontsize=7, rotation=0)
    ax.text((a + b) / 2, 0.22, "/".join(partners[:3]) + ("…" if len(partners) > 3 else ""), ha="center", fontsize=6, color="#2d3436")
ax.axvspan(TRUNC_END + 1, 1050, color="#c0392b", alpha=0.08)
ax.annotate("p.Leu737Ter\n(NMD likely; if translated, residues 737–1050 lost)", (737, 0.65), (600, 1.05), fontsize=8, ha="center",
            arrowprops=dict(arrowstyle="->", color="#c0392b"), color="#c0392b")
ax.annotate("p.Asn1002Lys\nAlphaMissense 0.92", (1002, 0.65), (960, 1.05), fontsize=8, ha="center",
            arrowprops=dict(arrowstyle="->", color="#d35400"), color="#d35400")
ax.set_xlim(0, 1060); ax.set_ylim(0, 1.3); ax.set_yticks([]); ax.set_xlabel("BUBR1 (BUB1B, UniProt O60566) residue")
ax.set_title("BUB1B protein map: functional regions, their interaction partners, and the proband's two variants", fontsize=10)
for s in ("top","right","left"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(R / "BUB1B_protein_map.png", dpi=170); plt.close(fig)

# ---------------------------------------------------------------- 6. Druggability: DGIdb + ChEMBL for altered proteins + partners
targets = sorted(set(altered) | set(diff.partner))
q = """query($n:[String!]!){ genes(names:$n){ nodes{ name interactions{ drug{ name approved } interactionScore interactionTypes{ type directionality } sources{ sourceDbName } } } } }"""
dg = S.post("https://dgidb.org/api/graphql", json={"query": q, "variables": {"n": targets}}, timeout=120).json()["data"]["genes"]["nodes"]
drows = []
for g in dg:
    for it in g["interactions"]:
        drows.append(dict(gene=g["name"], drug=it["drug"]["name"], approved=it["drug"]["approved"], score=round(it["interactionScore"], 3),
                          type=";".join(t["type"] for t in it["interactionTypes"]), sources=";".join(s["sourceDbName"] for s in it["sources"])))
dgi = pd.DataFrame(drows); dgi.to_csv(R / "dgidb_interactions.tsv", sep="\t", index=False)
sym2uni = dict(zip(hgnc.symbol, hgnc.uniprot))
crows = []
for g in targets:
    tj = get("https://www.ebi.ac.uk/chembl/api/data/target.json", params={"target_synonym__iexact": g, "organism": "Homo sapiens", "target_type": "SINGLE PROTEIN", "limit": 5}).json()
    tid = next((t["target_chembl_id"] for t in tj["targets"]), None)
    n_act = n_cpd = 0; max_phase = None
    if tid:
        a = get("https://www.ebi.ac.uk/chembl/api/data/activity.json", params={"target_chembl_id": tid, "pchembl_value__gte": 6, "limit": 1}).json()
        n_act = a["page_meta"]["total_count"]
        m = get("https://www.ebi.ac.uk/chembl/api/data/mechanism.json", params={"target_chembl_id": tid, "limit": 50}).json()["mechanisms"]
        n_cpd = len({x["molecule_chembl_id"] for x in m})
        if m:
            ph = [get("https://www.ebi.ac.uk/chembl/api/data/molecule/" + x["molecule_chembl_id"] + ".json").json().get("max_phase") for x in m[:10]]
            ph = [float(p) for p in ph if p is not None]; max_phase = max(ph) if ph else None
    sub = dgi[dgi.gene == g] if len(dgi) else dgi
    crows.append(dict(gene=g, altered_in_proband=g in altered, chembl_target=tid or "", chembl_potent_activities_pchembl_ge6=n_act,
                      chembl_mechanism_compounds=n_cpd, chembl_max_clinical_phase=max_phase,
                      dgidb_drugs=len(sub), dgidb_approved=int(sub.approved.sum()) if len(sub) else 0,
                      dgidb_top=";".join(sub.sort_values("score", ascending=False).drug.head(3)) if len(sub) else ""))
drug = pd.DataFrame(crows).sort_values(["altered_in_proband","chembl_potent_activities_pchembl_ge6"], ascending=[False, False])
drug.to_csv(R / "druggability.tsv", sep="\t", index=False)
print(drug.to_string(index=False))
print("Tier B done ->", R)
