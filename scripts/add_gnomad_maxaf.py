"""Append MAX_AF and MAX_AF_POPS columns to a VEP --tab file from a gnomAD lookup table (set,CHROM,POS,REF,ALT,AF,AF_grpmax,grpmax).
Usage: add_gnomad_maxaf.py vep_in.tsv gnomad_lookup.tsv vep_out.tsv"""
import sys, pandas as pd
vin, lk, vout = sys.argv[1:4]
g = pd.read_csv(lk, sep="\t", dtype=str)
for c in ("AF", "AF_grpmax"): g[c] = pd.to_numeric(g[c].replace(".", None), errors="coerce")
g["key"] = g.CHROM + "_" + g.POS + "_" + g.REF + "_" + g.ALT
g["best"] = g[["AF", "AF_grpmax"]].max(axis=1)
agg = g.groupby("key").apply(lambda d: pd.Series({"MAX_AF": d.best.max(),
      "MAX_AF_POPS": ";".join(f"{r.set}:{r.grpmax}" for r in d.itertuples() if pd.notna(r.best) and r.best == d.best.max())})).reset_index()
m = dict(zip(agg.key, zip(agg.MAX_AF, agg.MAX_AF_POPS)))
lines = open(vin).read().splitlines()
with open(vout, "w") as f:
    for l in lines:
        if l.startswith("##"): f.write(l + "\n"); continue
        if l.startswith("#Uploaded_variation"): f.write(l + "\tMAX_AF\tMAX_AF_POPS\n"); continue
        k = l.split("\t", 1)[0]; v = m.get(k)
        f.write(l + ("\t%.6g\t%s" % v if v and pd.notna(v[0]) else "\t-\t-") + "\n")
print(f"{len(m)} variants matched in gnomAD; wrote {vout}")
