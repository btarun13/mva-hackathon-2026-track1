"""Tier A analysis: genome-wide (coding) rare protein-altering variants -> gene burden -> HGNC -> enrichment (ORA with
coverage background + pre-ranked GSEA) -> STRING overlay on the checkpoint module -> DGIdb druggability.
Reads work/genome/vep_coding.tsv + coding.gt.tsv. Writes results/tierA_*.  Gene-level outputs only (no genotypes) are
safe to commit; the per-variant table stays in work/.
Run: uv run --no-project --with pandas,requests,gseapy,networkx,matplotlib python scripts/tier_a_analysis.py
"""
import io, sys, time
from pathlib import Path
import pandas as pd, requests, gseapy as gp
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from mva_track1 import read_vep_tab, consequence_class, clinvar_score, variant_score, zygosity, CORE_MVA_GENES
R = ROOT / "results"; W = ROOT / "work/genome"
S = requests.Session(); S.headers["User-Agent"] = "mva-hackathon-tier-a"
RARE_AF = 0.001

# ------------------------------------------------------------------ 1. load, QC, rarity, consequence
vep = read_vep_tab(str(W / "vep_coding.tsv"))
gt = pd.read_csv(W / "coding.gt.tsv", sep="\t", header=None,
                 names=["ID","CHROM","POS","REF","ALT","FILTER","GT","DP","GQ","AD"], dtype=str)
df = vep.merge(gt, left_on="Uploaded_variation", right_on="ID", how="inner")
n0 = len(df)
df = df[df["FILTER"].isin(["PASS", "."])]
df["DP"] = pd.to_numeric(df["DP"], errors="coerce"); df["GQ"] = pd.to_numeric(df["GQ"], errors="coerce")
df = df[(df["DP"].fillna(99) >= 10) & (df["GQ"].fillna(99) >= 20)]
df["zyg"] = df["GT"].map(zygosity); df = df[df["zyg"].isin(["het","hom"])]
df = df[df["BIOTYPE"].fillna("") == "protein_coding"]
background = sorted(df["SYMBOL"].dropna().unique())           # coverage-aware background: genes with any QC-pass coding variant
df["max_af"] = pd.to_numeric(df["MAX_AF"], errors="coerce")
df["cons_class"] = df["Consequence"].map(consequence_class)
df["clinvar_benign"] = df["ClinVar_CLNSIG"].map(lambda s: clinvar_score(s) is None)
prot = df[df["IMPACT"].isin(["HIGH","MODERATE"])]
rare = prot[(prot["max_af"].fillna(0) < RARE_AF) & ~prot["clinvar_benign"]].copy()
rare["vscore"] = [variant_score(r, 0.0) for _, r in rare.iterrows()]   # severity+rarity+ClinVar only, no gene prior
rare = rare.dropna(subset=["vscore"])
print(f"variants: {n0} annotated -> {len(df)} QC-pass protein-coding -> {len(prot)} HIGH/MODERATE -> {len(rare)} rare(<{RARE_AF}) non-benign")
print(f"background genes: {len(background)}")

# ------------------------------------------------------------------ 2. gene burden table
g = rare.groupby("SYMBOL").agg(n_rare=("ID","size"), n_hom=("zyg", lambda z: (z=="hom").sum()), n_het=("zyg", lambda z: (z=="het").sum()),
                               n_lof=("cons_class", lambda c: (c=="lof").sum()), max_vscore=("vscore","max"),
                               clinvar_plp=("ClinVar_CLNSIG", lambda s: int(s.fillna("").str.contains("athogenic").sum())),
                               ensg=("Gene","first"), top_hgvsp=("HGVSp", lambda s: s.dropna().iloc[0] if s.notna().any() else "")).reset_index()
g["recessive_compatible"] = (g.n_hom >= 1) | (g.n_het >= 2)
g["core_mva"] = g.SYMBOL.isin(CORE_MVA_GENES)
g = g.sort_values(["recessive_compatible","max_vscore","n_rare"], ascending=[False, False, False])
g.to_csv(R / "tierA_gene_burden.tsv", sep="\t", index=False)
rare_out = rare[["ID","SYMBOL","Gene","Consequence","IMPACT","HGVSc","HGVSp","max_af","MAX_AF_POPS","ClinVar_CLNSIG","am_class","am_pathogenicity","zyg","DP","GQ","vscore"]]
rare_out.to_csv(W / "rare_protein_altering_variants.tsv", sep="\t", index=False)   # patient-level: stays in work/
hits = sorted(g.SYMBOL); rec = sorted(g.loc[g.recessive_compatible, "SYMBOL"])
print(f"genes with >=1 rare protein-altering variant: {len(hits)} | recessive-compatible: {len(rec)} -> {rec[:40]}")

# ------------------------------------------------------------------ 3. HGNC for any Ensembl ID without a symbol
missing = rare.loc[rare["SYMBOL"].isna(), "Gene"].dropna().unique()
hg = []
for ensg in missing:
    d = S.get(f"https://rest.genenames.org/fetch/ensembl_gene_id/{ensg}", headers={"Accept":"application/json"}, timeout=30).json()["response"]["docs"]
    hg.append(dict(ensembl_gene_id=ensg, symbol=d[0]["symbol"] if d else "", name=d[0]["name"] if d else "(not in HGNC)"))
pd.DataFrame(hg).to_csv(R / "tierA_hgnc_unnamed.tsv", sep="\t", index=False)
print("Ensembl IDs without a VEP symbol:", len(missing))

# ------------------------------------------------------------------ 4. enrichment: ORA (hypergeometric, coverage background) + pre-ranked GSEA
LIBS = ["GO_Biological_Process_2023", "Reactome_2022", "KEGG_2021_Human", "Human_Phenotype_Ontology"]
libs = {l: gp.get_library(name=l, organism="Human") for l in LIBS}
ora = []
for label, genes in [("all_rare_hits", hits), ("recessive_compatible", rec)]:
    if len(genes) < 3: continue
    for l, gs in libs.items():
        e = gp.enrich(gene_list=genes, gene_sets=gs, background=background, outdir=None).results
        e.insert(0, "gene_list", label); e.insert(1, "library", l); ora.append(e)
ora = pd.concat(ora)
ora = ora[["gene_list","library","Term","Overlap","P-value","Adjusted P-value","Odds Ratio","Genes"]].sort_values(["gene_list","Adjusted P-value"])
ora.to_csv(R / "tierA_enrichment_ora.tsv", sep="\t", index=False)
sig = ora[ora["Adjusted P-value"] < 0.05]
print("ORA terms FDR<0.05:", sig.groupby("gene_list").size().to_dict())

rnk = pd.Series(0.0, index=background); rnk.update(g.set_index("SYMBOL")["max_vscore"]); rnk = rnk.sort_values(ascending=False)
rnk.to_csv(R / "tierA_gene_rank.rnk", sep="\t", header=False)
# Pre-ranked GSEA is not run: >95% of background genes tie at 0, so the ranking is arbitrary (see report).

# dot plot of top ORA terms per list
fig, axes = plt.subplots(1, 2, figsize=(16, 7))
for ax, label in zip(axes, ["all_rare_hits", "recessive_compatible"]):
    t = ora[ora.gene_list == label].sort_values("Adjusted P-value").head(15).iloc[::-1]
    if len(t):
        ax.barh(t.Term.str.slice(0, 60) + " [" + t.library.str.split("_").str[0] + "]", -t["Adjusted P-value"].map(lambda x: __import__("math").log10(max(x, 1e-300))), color="#2980b9")
        ax.axvline(-__import__("math").log10(0.05), color="#c0392b", ls="--", lw=1)
    ax.set_xlabel("-log10 FDR"); ax.set_title(f"ORA: {label} (n={len(hits) if label=='all_rare_hits' else len(rec)} genes)", fontsize=10); ax.tick_params(labelsize=7)
fig.tight_layout(); fig.savefig(R / "tierA_enrichment_ora.png", dpi=150); plt.close(fig)

# ------------------------------------------------------------------ 5. STRING: do recessive-compatible genes touch the checkpoint module?
SEED = ["BUB1B","BUB1","BUB3","MAD1L1","MAD2L1","CDC20","TRIP13","MAD2L1BP","CEP57","CENPE","KIF18A","CENATAC","TTK","PLK1","AURKB","KNL1","PPP2R5B","PPP2R5C","PPP2R5E","PPP2CA","CDC27","ANAPC1","CDK1","ZW10"]
query = sorted(set(SEED) | set(rec))
def string_tsv(endpoint, ids, **extra):
    r = S.post(f"https://string-db.org/api/tsv/{endpoint}", data={"identifiers": "\r".join(ids), "species": 9606, **extra}, timeout=300)
    t = pd.read_csv(io.StringIO(r.text), sep="\t")
    if "stringId" not in t.columns and "preferredName_A" not in t.columns:
        raise RuntimeError(f"STRING {endpoint} returned: {r.text[:300]}")
    return t
ids = string_tsv("get_string_ids", query, limit=1, echo_query=1)
sid = ids.drop_duplicates("queryItem")
net = string_tsv("network", list(sid.stringId), required_score=700, network_type="functional")
net = net.drop_duplicates(subset=["preferredName_A","preferredName_B"])
in_rec = lambda s: s.isin(rec) & ~s.isin(SEED)
link = net[(in_rec(net.preferredName_A) & net.preferredName_B.isin(SEED)) | (in_rec(net.preferredName_B) & net.preferredName_A.isin(SEED))]
link[["preferredName_A","preferredName_B","score","escore","dscore","tscore"]].to_csv(R / "tierA_string_links_to_checkpoint_module.tsv", sep="\t", index=False)
touching = sorted((set(link.preferredName_A) | set(link.preferredName_B)) - set(SEED))
print("recessive-compatible genes (outside the seed module) with a STRING edge >=0.7 into the checkpoint module:", touching)
print("BUB1B recessive-compatible genome-wide:", "BUB1B" in rec)

# ------------------------------------------------------------------ 6. DGIdb for recessive-compatible genes
q = """query($n:[String!]!){ genes(names:$n){ nodes{ name interactions{ drug{ name approved } interactionScore interactionTypes{ type } } } } }"""
drows = []
for i in range(0, len(rec), 50):
    nodes = S.post("https://dgidb.org/api/graphql", json={"query": q, "variables": {"n": rec[i:i+50]}}, timeout=120).json()["data"]["genes"]["nodes"]
    for gn in nodes:
        for it in gn["interactions"]:
            drows.append(dict(gene=gn["name"], drug=it["drug"]["name"], approved=it["drug"]["approved"], score=round(it["interactionScore"], 3), type=";".join(t["type"] for t in it["interactionTypes"])))
dgi = pd.DataFrame(drows); dgi.to_csv(R / "tierA_dgidb_recessive_genes.tsv", sep="\t", index=False)
print("DGIdb: recessive-compatible genes with any drug interaction:", dgi.gene.nunique() if len(dgi) else 0)
print("Tier A analysis done ->", R)
