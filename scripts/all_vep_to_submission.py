"""Write EVERY variant in work/vep.tsv as a row in the Track 1 submission CSV format (no filtering, no row cap).
Usage (from repo root): uv run --no-project --with pandas python scripts/all_vep_to_submission.py [USER_TAG]
Output: submissions2/<USER_TAG>_vep-rules-v1_all-variants.csv
"""
import csv, sys
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mva_track1 import (read_vep_tab, read_gene_list, consequence_class, rarity_score, clinvar_score,
                        spliceai_max, zygosity, CORE_MVA_GENES, PROBAND_ID, _s)

user_tag = sys.argv[1] if len(sys.argv) > 1 else "Tarun2712"
vep = read_vep_tab("work/vep.tsv")
gt = pd.read_csv("work/subset.gt.tsv", sep="\t", header=None,
                 names=["ID", "CHROM", "POS", "REF", "ALT", "FILTER", "GT", "DP", "GQ", "AD"], dtype=str)
df = vep.merge(gt, left_on="Uploaded_variation", right_on="ID", how="left")
listed = set(read_gene_list("genes.txt"))
df["max_af"] = pd.to_numeric(df["MAX_AF"], errors="coerce")
df["cons_class"] = df["Consequence"].map(consequence_class)
df["zyg"] = df["GT"].fillna("./.").map(zygosity)

def prior(sym):
    s = str(sym).upper()
    return 1.0 if s in CORE_MVA_GENES else (0.5 if s in listed else 0.1)

def score_all(r):
    """Same weights as mva_track1.variant_score but never drops a variant; returns (score, drop_reasons)."""
    reasons = []
    cv = clinvar_score(r.get("ClinVar_CLNSIG"))
    if cv is None:
        cv = 0.0; reasons.append("ClinVar benign")
    sai = spliceai_max(r.get("SpliceAI_pred"))
    am = r.get("am_pathogenicity")
    am = float(am) if am is not None and not pd.isna(am) else None
    cls = r["cons_class"]
    if cls == "lof": sev = 1.0
    elif cls == "missense": sev = 0.55 if am is None else 0.3 + 0.6 * am
    elif cls == "splice_region": sev = 0.2 + 0.7 * sai
    else:
        sev = 0.6 * sai
        if sev < 0.3: reasons.append("non-coding/synonymous without SpliceAI support")
    sev = max(sev, 0.9 * sai)
    rar = rarity_score(r["max_af"])
    if rar == 0.0: reasons.append("maxAF>1%")
    if r["zyg"] not in ("het", "hom"): reasons.append(f"GT={r.get('GT')}")
    if str(r.get("FILTER")) not in ("PASS", ".", "nan"): reasons.append(f"FILTER={r.get('FILTER')}")
    return 0.35 * sev + 0.20 * rar + 0.20 * cv + 0.25 * prior(r["SYMBOL"]), reasons

df["score"], df["drop"] = zip(*[score_all(r) for _, r in df.iterrows()])
df = df.sort_values(["score", "CHROM", "POS"], ascending=[False, True, True]).reset_index(drop=True)

top = df["score"].max()
rows = []
for i, v in df.iterrows():
    e = min(0.99, max(0.01, v["score"] / top * 0.95)) - i * 0.001
    note = (f"{_s(v,'SYMBOL')} {v['zyg']} {v['cons_class']} ({_s(v,'Consequence')}); "
            f"{_s(v,'HGVSc')} {_s(v,'HGVSp')}; ClinVar={_s(v,'ClinVar_CLNSIG')}; "
            f"maxAF={v['max_af'] if not pd.isna(v['max_af']) else 'absent'}")
    if v["drop"]:
        note += "; pipeline would drop: " + " & ".join(v["drop"])
    rows.append({
        "proband_id": PROBAND_ID, "chrom_1": v["CHROM"], "pos_1": int(v["POS"]), "ref_1": v["REF"], "alt_1": v["ALT"],
        "chrom_2": "", "pos_2": "", "ref_2": "", "alt_2": "",
        "epcr": f"{max(round(e, 3), 0.001):.3f}", "finding_type": "primary",
        "notes": note.replace(",", ";")[:250]})

cols = ["proband_id","chrom_1","pos_1","ref_1","alt_1","chrom_2","pos_2","ref_2","alt_2","epcr","finding_type","notes"]
out = Path(f"submissions2/{user_tag}_vep-rules-v1_all-variants.csv")
out.parent.mkdir(exist_ok=True)
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)
print(f"wrote {out} ({len(rows)} rows)")
