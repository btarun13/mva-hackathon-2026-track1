#!/usr/bin/env python3
"""
MVA Hackathon 2026 - Track 1 pipeline (VCF-first).

Subcommands (run in order):
  genes     -> build GRCh38 BED from a gene-symbol list (Ensembl REST)
  subset    -> bcftools norm + subset VCF to the BED, set IDs, extract genotypes
  annotate  -> run VEP (tab output) with gnomAD max AF, ClinVar, AlphaMissense, SpliceAI
  rank      -> filter, score, compound-het grouping, write submission CSV
  validate  -> parse the CSV exactly as the official scorer does

Example:
  python mva_track1.py genes    --genes genes.txt --vcf data/WGS_EX2312012_HGWCNDSX7.vcf.gz --out work/genes.bed
  python mva_track1.py subset   --vcf data/WGS_EX2312012_HGWCNDSX7.vcf.gz --bed work/genes.bed --ref GRCh38.fa --out work/subset.vcf.gz
  python mva_track1.py annotate --vcf work/subset.vcf.gz --out work/vep.tsv --clinvar clinvar_GRCh38.vcf.gz \
                                --am AlphaMissense_hg38.tsv.gz --spliceai-snv spliceai_snv.vcf.gz --spliceai-indel spliceai_indel.vcf.gz
  python mva_track1.py rank     --vep work/vep.tsv --gt work/subset.gt.tsv --genes genes.txt --out submissions/USERNAME_vep-rules-v1.csv
  python mva_track1.py validate --csv submissions/USERNAME_vep-rules-v1.csv
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

import pandas as pd

PROBAND_ID = "PROBAND01"  # the only accepted proband_id (see submit_track1.py)

# Core MVA genes (OMIM). Anything else in genes.txt is treated as a weaker DisGeNET-style association.
CORE_MVA_GENES = {"BUB1B", "CEP57", "TRIP13"}

LOF_TERMS = {
    "stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant",
    "start_lost", "stop_lost", "transcript_ablation",
}
MISSENSE_TERMS = {"missense_variant", "inframe_insertion", "inframe_deletion", "protein_altering_variant"}
SPLICE_REGION_TERMS = {"splice_region_variant", "splice_donor_5th_base_variant",
                       "splice_donor_region_variant", "splice_polypyrimidine_tract_variant"}


# ----------------------------------------------------------------------------- helpers
def run(cmd: list[str] | str, **kw) -> None:
    print("+", cmd if isinstance(cmd, str) else " ".join(cmd), file=sys.stderr)
    subprocess.run(cmd, shell=isinstance(cmd, str), check=True, **kw)


def need(tool: str) -> None:
    if shutil.which(tool) is None:
        sys.exit(f"'{tool}' not found on PATH")


def vcf_uses_chr_prefix(vcf: str) -> bool:
    out = subprocess.run(["bcftools", "view", "-h", vcf], capture_output=True, text=True, check=True).stdout
    contigs = [l for l in out.splitlines() if l.startswith("##contig")]
    return any("ID=chr" in l for l in contigs)


def read_gene_list(path: str) -> list[str]:
    return [l.strip().upper() for l in Path(path).read_text().splitlines() if l.strip() and not l.startswith("#")]


# ----------------------------------------------------------------------------- genes -> BED
def cmd_genes(a: argparse.Namespace) -> None:
    """Step 1/2: gene symbols (from OMIM/DisGeNET search) -> GRCh38 BED with padding."""
    genes = read_gene_list(a.genes)
    prefix = "chr" if vcf_uses_chr_prefix(a.vcf) else ""
    rows = []
    for g in genes:
        url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{g}?content-type=application/json;expand=0"
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                d = json.load(r)
        except Exception as e:  # noqa: BLE001
            print(f"  ! {g}: Ensembl lookup failed ({e}); add coordinates manually", file=sys.stderr)
            continue
        start = max(0, d["start"] - 1 - a.pad)  # BED is 0-based half-open
        end = d["end"] + a.pad
        rows.append((f"{prefix}{d['seq_region_name']}", start, end, g))
        print(f"  {g}\t{prefix}{d['seq_region_name']}:{d['start']}-{d['end']}", file=sys.stderr)
    rows.sort()
    with open(a.out, "w") as f:
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")
    print(f"wrote {a.out} ({len(rows)} genes)", file=sys.stderr)


# ----------------------------------------------------------------------------- subset VCF
def cmd_subset(a: argparse.Namespace) -> None:
    """Step 2: normalise, subset to BED, set ID=CHROM_POS_REF_ALT, dump genotypes."""
    need("bcftools"); need("tabix")
    out = a.out
    gt_out = out.replace(".vcf.gz", ".gt.tsv")
    norm = ["bcftools", "norm", "-m", "-both"] + (["-f", a.ref] if a.ref else []) + [a.vcf, "-Ou"]
    cmd = (
        " ".join(norm)
        + f" | bcftools view -T {a.bed} -Ou"
        + " | bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%ALT' -Oz -o " + out
    )
    run(cmd)
    run(["tabix", "-f", out])
    # genotype table joined later on ID
    run(
        "bcftools query -f '%ID\\t%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER[\\t%GT\\t%DP\\t%GQ\\t%AD]\\n' "
        + out + " > " + gt_out
    )
    n = sum(1 for _ in open(gt_out))
    print(f"wrote {out} and {gt_out} ({n} variants)", file=sys.stderr)


# ----------------------------------------------------------------------------- VEP
def cmd_annotate(a: argparse.Namespace) -> None:
    """Step 3: VEP tab output. Requires a VEP install/cache (or run inside the ensemblorg/ensembl-vep container)."""
    need(a.vep)
    cmd = [
        a.vep, "-i", a.vcf, "-o", a.out, "--tab", "--force_overwrite",
        "--assembly", "GRCh38",
        "--pick", "--symbol", "--canonical", "--hgvs", "--numbers",
    ]
    if a.gtf:  # cache-free mode: transcripts from a GTF, gnomAD from sliced custom VCFs
        cmd += ["--gtf", a.gtf]
    else:
        cmd += ["--cache", "--offline", "--dir_cache", a.cache, "--max_af", "--af_gnomade", "--af_gnomadg",
                "--fields", ("Uploaded_variation,Location,Allele,SYMBOL,Gene,Feature,CANONICAL,Consequence,IMPACT,"
                             "HGVSc,HGVSp,MAX_AF,MAX_AF_POPS,gnomADe_AF,gnomADg_AF,"
                             "ClinVar_CLNSIG,ClinVar_CLNREVSTAT,am_class,am_pathogenicity,SpliceAI_pred")]
    if a.fasta:
        cmd += ["--fasta", a.fasta]
    if a.clinvar:
        cmd += ["--custom", f"file={a.clinvar},short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT"]
    for name, path in (("gnomADe", a.gnomade), ("gnomADg", a.gnomadg)):
        if path:
            cmd += ["--custom", f"file={path},short_name={name},format=vcf,type=exact,coords=0,"
                                "fields=AF%AF_grpmax%nhomalt"]
    if a.am:
        cmd += ["--dir_plugins", a.plugins, "--plugin", f"AlphaMissense,file={a.am}"]
    if a.spliceai_snv and a.spliceai_indel:
        cmd += ["--plugin", f"SpliceAI,snv={a.spliceai_snv},indel={a.spliceai_indel}"]
    run(cmd)
    if a.gtf:
        add_max_af(a.out)
    print(f"wrote {a.out}", file=sys.stderr)


def add_max_af(path: str) -> None:
    """Cache-free mode has no VEP MAX_AF: use max of gnomAD v4.1 exome/genome AF_grpmax and AF."""
    lines = Path(path).read_text().splitlines()
    i = next(k for k, l in enumerate(lines) if l.startswith("#Uploaded_variation"))
    header = lines[i].split("\t")
    src = [c for c in ("gnomADe_AF_grpmax", "gnomADg_AF_grpmax", "gnomADe_AF", "gnomADg_AF") if c in header]
    out = lines[:i] + [lines[i] + "\tMAX_AF"]
    for l in lines[i + 1:]:
        f = l.split("\t")
        afs = [float(f[header.index(c)]) for c in src if f[header.index(c)] not in ("-", "", ".")]
        out.append(l + "\t" + (f"{max(afs):.6g}" if afs else "-"))
    Path(path).write_text("\n".join(out) + "\n")


def read_vep_tab(path: str) -> pd.DataFrame:
    header = None
    with open(path) as f:
        for line in f:
            if line.startswith("#Uploaded_variation"):
                header = line[1:].rstrip("\n").split("\t")
                break
    if header is None:
        sys.exit("VEP header line not found")
    df = pd.read_csv(path, sep="\t", comment="#", names=header, na_values=["-"], dtype=str)
    return df


# ----------------------------------------------------------------------------- scoring
def consequence_class(cons: str) -> str:
    terms = set(str(cons).split(","))
    if terms & LOF_TERMS:
        return "lof"
    if terms & MISSENSE_TERMS:
        return "missense"
    if terms & SPLICE_REGION_TERMS:
        return "splice_region"
    return "other"


def rarity_score(max_af: float | None) -> float:
    if max_af is None or math.isnan(max_af):
        return 1.0  # absent from gnomAD
    if max_af < 1e-4:
        return 0.9
    if max_af < 1e-3:
        return 0.7
    if max_af < 1e-2:
        return 0.4
    return 0.0


def clinvar_score(sig: str | None) -> float | None:
    """Returns None for 'drop this variant' (benign), else a 0-1 contribution."""
    s = str(sig or "").lower()
    if not s or s == "nan":
        return 0.3  # unreported
    if "pathogenic" in s and "conflicting" not in s and "benign" not in s:
        return 1.0
    if "conflicting" in s or "uncertain" in s:
        return 0.4
    if "benign" in s:
        return None
    return 0.3


def spliceai_max(pred: str | None) -> float:
    # SpliceAI_pred format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
    if not pred or str(pred) == "nan":
        return 0.0
    best = 0.0
    for p in str(pred).split(","):
        parts = p.split("|")
        if len(parts) >= 6:
            try:
                best = max(best, *[float(x) for x in parts[2:6] if x not in ("", ".")])
            except ValueError:
                pass
    return best


def variant_score(row: pd.Series, gene_prior: float) -> float | None:
    cls = row["cons_class"]
    cv = clinvar_score(row.get("ClinVar_CLNSIG"))
    if cv is None:
        return None  # ClinVar benign/likely benign -> drop
    sai = spliceai_max(row.get("SpliceAI_pred"))
    am = row.get("am_pathogenicity")
    am = float(am) if am not in (None, "nan") and not pd.isna(am) else None

    if cls == "lof":
        sev = 1.0
    elif cls == "missense":
        sev = 0.55 if am is None else 0.3 + 0.6 * am  # AlphaMissense rescales 0.3-0.9
    elif cls == "splice_region":
        sev = 0.2 + 0.7 * sai
    else:
        sev = 0.6 * sai  # intronic/synonymous only survive if SpliceAI says splice-disrupting
        if sev < 0.3:
            return None
    sev = max(sev, 0.9 * sai)  # strong SpliceAI hit upgrades any class

    rar = rarity_score(row["max_af"])
    if rar == 0.0:
        return None  # > 1 % in some population: too common for a recessive Mendelian disease

    # weights: severity 0.35, rarity 0.20, ClinVar 0.20, gene prior 0.25
    return 0.35 * sev + 0.20 * rar + 0.20 * cv + 0.25 * gene_prior


def zygosity(gt: str) -> str:
    g = gt.replace("|", "/")
    if g in ("1/1",):
        return "hom"
    if g in ("0/1", "1/0", "1/2", "2/1"):
        return "het"
    return "other"


def cmd_rank(a: argparse.Namespace) -> None:
    """Steps 4-6 + extras A/B/C: filter, score, group compound hets, write submission."""
    vep = read_vep_tab(a.vep)
    gt = pd.read_csv(a.gt, sep="\t", header=None,
                     names=["ID", "CHROM", "POS", "REF", "ALT", "FILTER", "GT", "DP", "GQ", "AD"], dtype=str)
    df = vep.merge(gt, left_on="Uploaded_variation", right_on="ID", how="inner")
    print(f"{len(df)} annotated variants in candidate genes", file=sys.stderr)

    # --- basic QC
    df = df[df["FILTER"].isin(["PASS", "."])]
    df["DP"] = pd.to_numeric(df["DP"], errors="coerce")
    df["GQ"] = pd.to_numeric(df["GQ"], errors="coerce")
    df = df[(df["DP"].fillna(99) >= a.min_dp) & (df["GQ"].fillna(99) >= a.min_gq)]
    df["zyg"] = df["GT"].map(zygosity)
    df = df[df["zyg"].isin(["het", "hom"])]

    # --- gene prior (core MVA 1.0, other listed 0.5), optional Exomiser gene score boost
    listed = set(read_gene_list(a.genes)) if a.genes else set()
    exo = {}
    if a.exomiser_genes:
        eg = pd.read_csv(a.exomiser_genes, sep="\t")
        col = next((c for c in eg.columns if "COMBINED" in c.upper()), None)
        sym = next((c for c in eg.columns if "SYMBOL" in c.upper()), None)
        if col and sym:
            exo = dict(zip(eg[sym].str.upper(), eg[col].astype(float)))

    def prior(sym: str) -> float:
        s = str(sym).upper()
        base = 1.0 if s in CORE_MVA_GENES else (0.5 if s in listed else 0.1)
        return max(base, 0.5 * base + 0.5 * exo.get(s, 0.0)) if exo else base

    # --- rarity + consequence + score
    df["max_af"] = pd.to_numeric(df["MAX_AF"], errors="coerce")
    df["cons_class"] = df["Consequence"].map(consequence_class)
    df["gene_prior"] = df["SYMBOL"].map(prior)
    df["score"] = [variant_score(r, r["gene_prior"]) for _, r in df.iterrows()]
    df = df.dropna(subset=["score"]).sort_values("score", ascending=False)
    print(f"{len(df)} variants after rarity/consequence/ClinVar filters", file=sys.stderr)
    df.to_csv(a.table, sep="\t", index=False)
    print(f"wrote full scored table -> {a.table}", file=sys.stderr)

    # --- Extra A: group by gene -> rows (pair rows for recessive-compatible genes)
    rows: list[dict] = []
    for gene, g in df.groupby("SYMBOL", sort=False):
        g = g.sort_values("score", ascending=False)
        homs = g[g["zyg"] == "hom"]
        hets = g[g["zyg"] == "het"]
        for _, v in homs.iterrows():  # homozygous: single-variant row, full recessive bonus
            rows.append(_row(v, None, v["score"] + 0.15, f"{gene} homozygous {v['cons_class']}"))
        if len(hets) >= 2:  # compound-het candidates: best pairs on ONE row each
            for v1, v2 in itertools.combinations([r for _, r in hets.head(a.max_hets_per_gene).iterrows()], 2):
                pair_score = 0.5 * (v1["score"] + v2["score"]) + 0.15
                rows.append(_row(v1, v2, pair_score,
                                 f"{gene} compound-het candidate (phase unknown, single-sample VCF)"))
        for _, v in hets.iterrows():  # lone hets kept, but ranked below any recessive-compatible row
            rows.append(_row(v, None, v["score"] - 0.10, f"{gene} single het {v['cons_class']}"))

    rows.sort(key=lambda r: -r["_score"])
    rows = rows[: a.max_rows]

    # --- Extra B: map score -> EPCR in (0,1], strictly decreasing so thresholds are distinct
    if rows:
        top = rows[0]["_score"]
        for i, r in enumerate(rows):
            e = min(0.99, max(0.01, r["_score"] / top * 0.95))
            r["epcr"] = round(e - i * 0.001, 3)  # tiny monotone decrement keeps order unambiguous
            r["epcr"] = max(r["epcr"], 0.01)

    cols = ["proband_id", "chrom_1", "pos_1", "ref_1", "alt_1", "chrom_2", "pos_2", "ref_2", "alt_2",
            "epcr", "finding_type", "notes"]
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})
    print(f"wrote submission -> {a.out} ({len(rows)} rows)", file=sys.stderr)
    for r in rows:
        print(f"  {r['epcr']:.3f}  {r['chrom_1']}:{r['pos_1']} {r['ref_1']}>{r['alt_1']}"
              + (f"  +  {r['chrom_2']}:{r['pos_2']} {r['ref_2']}>{r['alt_2']}" if r['chrom_2'] else "")
              + f"   {r['notes']}", file=sys.stderr)


def _s(v: pd.Series, key: str) -> str:
    x = v.get(key)
    return "NA" if x is None or (isinstance(x, float) and math.isnan(x)) or str(x) == "nan" else str(x)


def _row(v1: pd.Series, v2: pd.Series | None, score: float, note: str) -> dict:
    r = {
        "proband_id": PROBAND_ID,
        "chrom_1": v1["CHROM"], "pos_1": int(v1["POS"]), "ref_1": v1["REF"], "alt_1": v1["ALT"],
        "chrom_2": "", "pos_2": "", "ref_2": "", "alt_2": "",
        "finding_type": "primary",
        "notes": f"{note}; {_s(v1, 'HGVSc')} {_s(v1, 'HGVSp')}; ClinVar={_s(v1, 'ClinVar_CLNSIG')}; "
                 f"maxAF={v1['max_af'] if not pd.isna(v1['max_af']) else 'absent'}",
        "_score": score,
    }
    if v2 is not None:
        r.update({"chrom_2": v2["CHROM"], "pos_2": int(v2["POS"]), "ref_2": v2["REF"], "alt_2": v2["ALT"]})
        r["notes"] += (f" | {_s(v2, 'HGVSc')} {_s(v2, 'HGVSp')}; ClinVar={_s(v2, 'ClinVar_CLNSIG')}; "
                       f"maxAF={v2['max_af'] if not pd.isna(v2['max_af']) else 'absent'}")
    r["notes"] = r["notes"].replace(",", ";")[:250]
    return r


# ----------------------------------------------------------------------------- validate (mirrors evaluation.py)
def cmd_validate(a: argparse.Namespace) -> None:
    with open(a.csv, newline="") as f:
        reader = csv.DictReader(f)
        required = {"proband_id", "chrom_1", "pos_1", "ref_1", "alt_1", "epcr"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            sys.exit(f"missing columns: {missing}")
        rows = list(reader)
    if not rows:
        sys.exit("no rows")
    if len(rows) > 10:
        sys.exit(f"{len(rows)} rows; max is 10")
    prev = None
    for i, r in enumerate(rows, 1):
        if r["proband_id"].strip() != PROBAND_ID:
            sys.exit(f"row {i}: proband_id must be {PROBAND_ID}")
        if not r["chrom_1"].startswith("chr"):
            print(f"row {i}: chrom_1='{r['chrom_1']}' lacks 'chr' prefix - scorer compares strings exactly", file=sys.stderr)
        int(r["pos_1"])
        if r.get("chrom_2"):
            int(r["pos_2"])
            if r["chrom_2"] == r["chrom_1"] and r["pos_2"] == r["pos_1"]:
                sys.exit(f"row {i}: pair contains the same variant twice")
        e = float(r["epcr"])
        if not 0 < e <= 1:
            sys.exit(f"row {i}: epcr {e} outside (0,1]")
        if prev is not None and e > prev:
            print(f"row {i}: epcr not sorted descending (scorer re-sorts, but pre-sort anyway)", file=sys.stderr)
        prev = e
        ft = (r.get("finding_type") or "primary").strip().lower()
        if ft not in ("primary", "secondary"):
            sys.exit(f"row {i}: finding_type must be primary/secondary")
    print(f"OK: {len(rows)} rows, parses like the official scorer", file=sys.stderr)


# ----------------------------------------------------------------------------- CLI
def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest="cmd", required=True)

    g = sp.add_parser("genes"); g.add_argument("--genes", required=True); g.add_argument("--vcf", required=True)
    g.add_argument("--out", required=True); g.add_argument("--pad", type=int, default=5000); g.set_defaults(fn=cmd_genes)

    s = sp.add_parser("subset"); s.add_argument("--vcf", required=True); s.add_argument("--bed", required=True)
    s.add_argument("--ref", help="reference FASTA for left-alignment (optional)"); s.add_argument("--out", required=True)
    s.set_defaults(fn=cmd_subset)

    v = sp.add_parser("annotate"); v.add_argument("--vcf", required=True); v.add_argument("--out", required=True)
    v.add_argument("--vep", default="vep"); v.add_argument("--cache", default=str(Path.home() / ".vep"))
    v.add_argument("--fasta"); v.add_argument("--clinvar"); v.add_argument("--am")
    v.add_argument("--spliceai-snv"); v.add_argument("--spliceai-indel")
    v.add_argument("--gtf", help="cache-free mode: bgzipped+tabixed GTF instead of the VEP cache")
    v.add_argument("--gnomade"); v.add_argument("--gnomadg"); v.add_argument("--plugins", default="/plugins")
    v.set_defaults(fn=cmd_annotate)

    r = sp.add_parser("rank"); r.add_argument("--vep", required=True); r.add_argument("--gt", required=True)
    r.add_argument("--genes"); r.add_argument("--exomiser-genes", help="Exomiser *.genes.tsv for phenotype prior")
    r.add_argument("--out", required=True); r.add_argument("--table", default="work/scored_variants.tsv")
    r.add_argument("--min-dp", type=int, default=10); r.add_argument("--min-gq", type=int, default=20)
    r.add_argument("--max-hets-per-gene", type=int, default=3); r.add_argument("--max-rows", type=int, default=10)
    r.set_defaults(fn=cmd_rank)

    c = sp.add_parser("validate"); c.add_argument("--csv", required=True); c.set_defaults(fn=cmd_validate)

    a = p.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()
