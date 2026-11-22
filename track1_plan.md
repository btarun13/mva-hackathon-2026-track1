# Track 1 – Variant Prediction Plan
**Rare Disease, Real Kid: MVA Hackathon 2026**

All code lives in one script, `mva_track1.py` (subcommands `genes`, `subset`, `annotate`, `rank`, `validate`). Each step below shows the command that runs it and the core logic.

## Inputs provided
- Child's VCF (GRCh38) — primary input: `WGS_EX2312012_HGWCNDSX7.vcf.gz` + `.tbi`
- Child's FASTQ (8 files, ~85 GB) — only used in Plan B
- Disease: mosaic variegated aneuploidy (MVA)
- Clinical symptom description: `Challenge_Clinical_Phenotype_1.docx`

## Scoring to keep in mind (from the Space's `evaluation.py`)
- **Rank points**: rank 1 → 100, ranks 2–3 → 50, ranks 4–5 → 25, ranks 6–10 → 10, else 0. Half credit if a row contains only one of the two true compound-het variants.
- **A compound-het pair must be on ONE row** (`chrom_1..alt_1` + `chrom_2..alt_2`). Two separate rows can never earn more than half credit. Pair order doesn't matter.
- Matching is an exact string match on `(chrom, pos, ref, alt)` — so use `chr15`, not `15`, and GRCh38 coordinates.
- **F-max**: computed at the individual-variant level over the distinct `epcr` thresholds. If the top row alone is the correct pair, F-max = 1.0 regardless of what sits below it. Give the top row a unique highest `epcr`.
- `proband_id` must be `PROBAND01`; max 10 rows; `epcr` in (0, 1]; `finding_type` is `primary`/`secondary`.
- **Methods write-up**: judged separately for scientific rigor. Submission limit: 6.

## Setup
```bash
pip install pandas huggingface_hub
# bcftools + tabix via conda/apt; VEP via the ensemblorg/ensembl-vep container or a local install with the GRCh38 cache
hf download SageBio/mva-hackathon-2026-data \
  WGS_EX2312012_HGWCNDSX7.vcf.gz WGS_EX2312012_HGWCNDSX7.vcf.gz.tbi \
  Challenge_Clinical_Phenotype_1.docx README.md --repo-type dataset --local-dir ./data
bcftools view -h data/WGS_EX2312012_HGWCNDSX7.vcf.gz | grep -E "##reference|##contig=<ID=chr1," | head   # confirm GRCh38 + chr prefix
```

---

## Core pipeline

### Step 1 – Disease gene list (OMIM + DisGeNET)
Search "mosaic variegated aneuploidy" in [OMIM](https://omim.org) and [DisGeNET](https://www.disgenet.com). Put one symbol per line in `genes.txt`. Core genes `BUB1B`, `TRIP13`, `CEP57` get gene-prior 1.0 in scoring; anything else in the file gets 0.5.

```text
# genes.txt
BUB1B
CEP57
TRIP13
CENPE
KIF18A
# + DisGeNET associations
```

### Step 2 – Subset the VCF to candidate genes
`genes` looks up GRCh38 coordinates from Ensembl REST and writes a padded BED (chr prefix chosen to match the VCF). `subset` normalises (split multiallelics, left-align), restricts to the BED, sets `ID = CHROM_POS_REF_ALT` so VEP output can be joined back, and dumps genotypes.

```bash
python mva_track1.py genes  --genes genes.txt --vcf data/WGS_EX2312012_HGWCNDSX7.vcf.gz --out work/genes.bed
python mva_track1.py subset --vcf data/WGS_EX2312012_HGWCNDSX7.vcf.gz --bed work/genes.bed --ref GRCh38.fa --out work/subset.vcf.gz
# -> work/subset.vcf.gz and work/subset.gt.tsv (ID CHROM POS REF ALT FILTER GT DP GQ AD)
```

Core of `subset`:
```python
cmd = (f"bcftools norm -m -both -f {ref} {vcf} -Ou"
       f" | bcftools view -R {bed} -Ou"
       " | bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%ALT' -Oz -o work/subset.vcf.gz")
run(cmd); run(["tabix", "-f", "work/subset.vcf.gz"])
run("bcftools query -f '%ID\\t%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER[\\t%GT\\t%DP\\t%GQ\\t%AD]\\n' "
    "work/subset.vcf.gz > work/subset.gt.tsv")
```

### Step 3 – Annotate with VEP
`annotate` runs VEP in tab mode with `--pick` (one consequence per variant), `--max_af` (max AF across gnomAD/1000G populations), gnomAD exome+genome AFs, ClinVar as a custom VCF, and (Extra D) AlphaMissense + SpliceAI plugins when the files are given.

```bash
python mva_track1.py annotate --vcf work/subset.vcf.gz --out work/vep.tsv \
  --cache ~/.vep --fasta GRCh38.fa \
  --clinvar clinvar_GRCh38.vcf.gz \
  --am AlphaMissense_hg38.tsv.gz \
  --spliceai-snv spliceai_scores.raw.snv.hg38.vcf.gz --spliceai-indel spliceai_scores.raw.indel.hg38.vcf.gz
```

VEP flags used:
```python
cmd = [vep, "-i", vcf, "-o", out, "--tab", "--force_overwrite",
       "--assembly", "GRCh38", "--cache", "--offline", "--dir_cache", cache,
       "--pick", "--symbol", "--canonical", "--hgvs", "--numbers",
       "--max_af", "--af_gnomade", "--af_gnomadg",
       "--fields", "Uploaded_variation,Location,Allele,SYMBOL,Gene,Feature,CANONICAL,Consequence,IMPACT,"
                   "HGVSc,HGVSp,MAX_AF,MAX_AF_POPS,gnomADe_AF,gnomADg_AF,"
                   "ClinVar_CLNSIG,ClinVar_CLNREVSTAT,am_class,am_pathogenicity,SpliceAI_pred",
       "--custom", f"{clinvar},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT",
       "--plugin", f"AlphaMissense,file={am}",
       "--plugin", f"SpliceAI,snv={spliceai_snv},indel={spliceai_indel}"]
```
Output per variant: gene, consequence, HGVS, ClinVar status, gnomAD max-population AF, AlphaMissense, SpliceAI.

### Step 4 – Rarity filter (gnomAD)
Uses `MAX_AF`. Max AF > 1 % → dropped (MVA is autosomal recessive). Absent from gnomAD scores highest.

```python
def rarity_score(max_af):
    if max_af is None or math.isnan(max_af): return 1.0   # absent from gnomAD
    if max_af < 1e-4: return 0.9
    if max_af < 1e-3: return 0.7
    if max_af < 1e-2: return 0.4
    return 0.0                                              # > 1 %: variant dropped
```

### Step 5 – Consequence filter
LoF and missense kept; synonymous/deep-intronic dropped **unless** SpliceAI says they disrupt splicing (that is the one case Extra D rescues). ClinVar benign/likely benign → dropped.

```python
LOF_TERMS      = {"stop_gained","frameshift_variant","splice_acceptor_variant","splice_donor_variant",
                  "start_lost","stop_lost","transcript_ablation"}
MISSENSE_TERMS = {"missense_variant","inframe_insertion","inframe_deletion","protein_altering_variant"}

def consequence_class(cons):
    t = set(cons.split(","))
    if t & LOF_TERMS: return "lof"
    if t & MISSENSE_TERMS: return "missense"
    if t & SPLICE_REGION_TERMS: return "splice_region"
    return "other"
```

Also applied: `FILTER` = PASS, `DP ≥ 10`, `GQ ≥ 20`, genotype het or hom.

### Step 6 – Rank and submit
`rank` does Steps 4–5, scoring (Extra B), compound-het grouping (Extra A), and writes the CSV in the exact scorer format. `validate` re-parses it the way `evaluation.py` does.

```bash
python mva_track1.py rank --vep work/vep.tsv --gt work/subset.gt.tsv --genes genes.txt \
  --out submissions/<hf-username>_vep-rules-v1.csv
python mva_track1.py validate --csv submissions/<hf-username>_vep-rules-v1.csv
```

Ordering follows the plan: ClinVar P/LP > LoF > missense, then rarity — expressed as a weighted score so AlphaMissense/SpliceAI and gene prior can break ties:
```python
def variant_score(row, gene_prior):
    cv = clinvar_score(row["ClinVar_CLNSIG"])        # P/LP 1.0, VUS/conflicting 0.4, none 0.3, benign -> None (drop)
    if cv is None: return None
    sai = spliceai_max(row["SpliceAI_pred"])         # max of DS_AG/AL/DG/DL
    am  = row["am_pathogenicity"]                    # AlphaMissense, may be missing
    cls = row["cons_class"]
    if   cls == "lof":           sev = 1.0
    elif cls == "missense":      sev = 0.55 if am is None else 0.3 + 0.6*am
    elif cls == "splice_region": sev = 0.2 + 0.7*sai
    else:                                            # synonymous / intronic
        sev = 0.6*sai
        if sev < 0.3: return None
    sev = max(sev, 0.9*sai)                          # strong SpliceAI hit upgrades any class
    rar = rarity_score(row["max_af"])
    if rar == 0.0: return None
    return 0.35*sev + 0.20*rar + 0.20*cv + 0.25*gene_prior
```

### Plan B (only if the score is poor)
Align FASTQs to GRCh38 and re-call, then rerun from Step 2 with the new VCF (`--vcf work/planB.vcf.gz`).
```bash
bwa-mem2 index GRCh38.fa
for L in 001 002 003 004; do
  bwa-mem2 mem -t 16 -R "@RG\tID:L$L\tSM:PROBAND01\tPL:ILLUMINA" GRCh38.fa \
    data/WGS_EX2312012_HGWCNDSX7_S16_L${L}_R1_001.fastq.gz data/WGS_EX2312012_HGWCNDSX7_S16_L${L}_R2_001.fastq.gz \
  | samtools sort -@8 -o work/L$L.bam
done
samtools merge -@8 work/proband.bam work/L00?.bam && samtools index work/proband.bam
docker run -v $PWD:/d google/deepvariant:1.6.1 /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS --ref=/d/GRCh38.fa --reads=/d/work/proband.bam --output_vcf=/d/work/planB.vcf.gz --num_shards=16
```

---

## Extras required by the scoring rules

### A. Compound-heterozygote grouping (in `rank`)
Surviving variants are grouped by gene. Homozygous hits get a recessive bonus. Any two heterozygous hits in the same gene are emitted together **on one row** (chrom_1/chrom_2), scored as their mean + bonus; lone hets are kept but pushed below any recessive-compatible row. Phase is unknown from a single-sample VCF — the notes say so.

```python
for gene, g in df.groupby("SYMBOL"):
    g = g.sort_values("score", ascending=False)
    homs, hets = g[g.zyg == "hom"], g[g.zyg == "het"]
    for _, v in homs.iterrows():
        rows.append(_row(v, None, v.score + 0.15, f"{gene} homozygous"))
    if len(hets) >= 2:
        for v1, v2 in itertools.combinations(list(hets.head(3).itertuples()), 2):
            rows.append(_row(v1, v2, 0.5*(v1.score + v2.score) + 0.15,
                             f"{gene} compound-het candidate (phase unknown)"))
    for _, v in hets.iterrows():
        rows.append(_row(v, None, v.score - 0.10, f"{gene} single het"))
rows.sort(key=lambda r: -r["_score"]); rows = rows[:10]
```

### B. Confidence score and cut-off (in `rank`)
Row scores are mapped to `epcr` in (0, 1], strictly decreasing so every row is its own F-max threshold. Because F-max is a maximum over thresholds, extra lower-ranked rows cost nothing — the thing that matters is that the top row is right and unique. `--max-rows` (default 10) caps the list.

```python
top = rows[0]["_score"]
for i, r in enumerate(rows):
    r["epcr"] = max(0.01, round(min(0.99, r["_score"]/top*0.95) - i*0.001, 3))
```

### C. Use the symptom description
Convert the docx symptoms to HPO terms, run Exomiser genome-wide, and feed its gene-level combined score into the gene prior:
```bash
# exomiser-config.yml: hpoIds from the docx, vcf: data/WGS_...vcf.gz, genomeAssembly: hg38,
#   inheritanceModes AUTOSOMAL_RECESSIVE + AUTOSOMAL_DOMINANT, frequency filter 1%, outputFormats [TSV_GENE, TSV_VARIANT, HTML]
java -jar exomiser-cli.jar --analysis exomiser-config.yml
python mva_track1.py rank ... --exomiser-genes results/proband.genes.tsv --out submissions/<user>_vep-rules-exomiser-v2.csv
```
```python
def prior(sym):                      # inside rank
    base = 1.0 if sym in CORE_MVA_GENES else (0.5 if sym in listed else 0.1)
    return max(base, 0.5*base + 0.5*exomiser_combined.get(sym, 0.0)) if exomiser_combined else base
```
Exomiser's own top-10 is also a sanity check: it should agree with the manual panel result.

### D. Optional deep-learning annotations
Already wired in Step 3 via `--am` and `--spliceai-*`. AlphaMissense rescales missense severity (0.3–0.9); SpliceAI ≥ 0.5 rescues otherwise-filtered synonymous/intronic variants and upgrades splice-region ones.

---

## Submission strategy
1. Run Steps 1–6 without `--exomiser-genes` → `..._vep-rules-v1.csv`. Submit, read rank points / F-max / match type (full vs partial tells you whether the pair is on one row).
2. If partial match: the second allele is in the filtered-out set — relax `--min-dp`/`--min-gq`, raise `--max-hets-per-gene`, or check for an SV/CNV in the same gene. Resubmit.
3. Add C (Exomiser prior) and D → `..._v2.csv`. Resubmit.
4. Keep at least two attempts in reserve for final tuning and the write-up's final list.

Each submission needs: the CSV (username + method in filename), a Markdown/PDF report, and a GitHub URL. If an LLM was used, add one line stating provider/tier/data-handling setting.

## Suggested tool stack
- `bcftools`, `tabix` – VCF normalisation, subsetting, genotype extraction
- VEP (GRCh38 cache) with ClinVar custom track, AlphaMissense and SpliceAI plugins
- Exomiser – phenotype-driven gene prior
- Python / pandas – `mva_track1.py`
- `bwa-mem2`, `samtools`, DeepVariant – Plan B only
