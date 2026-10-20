# MVA Hackathon 2026 — Track 1 report

**Team / user:** Tarun2712
**Model 1:** `vep-rules-v1` (submission file `Tarun2712_vep-rules-v1.csv`)
**Code:** https://github.com/btarun13/mva-hackathon-2026-track1 (`scripts/run_pipeline.sh` reproduces the submission byte-for-byte)

## Result

| Rank | EPCR | Variant(s) (GRCh38) | Gene | Interpretation |
|---|---|---|---|---|
| 1 | 0.950 | chr15:40209701 T>G **+** chr15:40220612 T>G | *BUB1B* | Compound-heterozygous candidate: c.2210T>G p.(Leu737Ter) + c.3006T>G p.(Asn1002Lys) |
| 2 | 0.807 | chr15:40209701 T>G | *BUB1B* | p.(Leu737Ter) alone (fallback) |
| 3 | 0.631 | chr15:40220612 T>G | *BUB1B* | p.(Asn1002Lys) alone (fallback) |

Transcript ENST00000287598 (canonical *BUB1B*).

| Evidence | c.2210T>G p.(Leu737Ter) | c.3006T>G p.(Asn1002Lys) |
|---|---|---|
| Consequence | stop_gained (truncates before the C-terminal kinase domain) | missense, C-terminal kinase domain |
| ClinVar (2026-10-19) | Pathogenic/Likely pathogenic | not reported |
| gnomAD v4.1 | exome AF 7.9e-5, genome AF 3.3e-5 | exome AF 6.8e-7, absent from genomes |
| AlphaMissense | — | 0.923 (likely_pathogenic) |
| Genotype (DP, AD ref/alt, GQ) | 0/1 (46, 21/25, 99), PASS | 0/1 (28, 15/13, 99), PASS |

Biallelic *BUB1B* variants cause mosaic variegated aneuploidy type 1 (MVA1, OMIM 257300). The type of cancer recorded for the proband (HPO HP:0002859) and the pre- and postnatal growth restriction are typical of *BUB1B*-associated MVA, more so than of *CEP57* (MVA2) or *TRIP13* (MVA3). This agrees with the purely genotype-driven ranking. The phenotype was **not** used to score variants (see Limitations).

## Approach

The approach is a candidate-gene panel with rules-based scoring. Recessive inheritance is explicitly modelled, and compound-het pairs are emitted together on one row. The pipeline is the single script `mva_track1.py`, run by `scripts/run_pipeline.sh`.

1. **Gene panel (`genes.txt`).** Core MVA genes from OMIM (*BUB1B*, *CEP57*, *TRIP13*) get gene prior 1.0. Related mitotic or spindle-assembly-checkpoint genes (*CENPE*, *KIF18A*, *CENATAC*, *BUB1*) get gene prior 0.5. Regions are the gene bodies from Ensembl 113 ±5 kb, so that promoter and UTR variation is included.
2. **VCF preparation.** The patient VCF (GATK, bare contig names) is restricted to the panel regions and renamed to `chr*` contigs. It is then normalised with `bcftools norm -m -both -f` (multiallelic records split, indels left-aligned). IDs are set to `CHROM_POS_REF_ALT` for a lossless join. The result is 252 records.
3. **Annotation.** VEP 113 (Docker, `--pick --canonical --hgvs`) runs in cache-free mode. It uses:
   - the Ensembl 113 GTF and FASTA;
   - ClinVar through `--custom`, exact allele match;
   - gnomAD v4.1 exome and genome sites through `--custom`, exact match on AF, AF_grpmax and nhomalt;
   - the AlphaMissense plugin.

   `MAX_AF` = max(AF_grpmax, AF) over the gnomAD exome and genome sets. All annotation sources were sliced to the panel regions by tabix region queries, about 250 MB in total.
4. **Filters.**
   - QC: FILTER = PASS, DP ≥ 10, GQ ≥ 20, het or hom-alt.
   - Rarity: MAX_AF < 1%; absent from gnomAD counts as rarest.
   - ClinVar: (likely) benign variants are removed.
   - Consequence: LoF, missense/in-frame and splice-region variants are kept. Synonymous, intronic and UTR variants survive only with SpliceAI ≥ 0.5, and SpliceAI was not available in this run.
   - Counts at each stage: 252 → 249 (PASS) → 244 (DP/GQ) → 10 (rare) → 2 (coding/splice).
5. **Variant score.** Score = 0.35·severity + 0.20·rarity + 0.20·ClinVar + 0.25·gene prior.
   - Severity: LoF 1.0; missense 0.3 + 0.6·AlphaMissense (0.55 if no AlphaMissense score).
   - Rarity: 1.0 if absent from gnomAD, 0.9 below 1e-4, 0.7 below 1e-3, 0.4 below 1e-2.
   - ClinVar: P/LP 1.0, VUS or conflicting 0.4, unreported 0.3.
6. **Recessive grouping.** Within each gene:
   - homozygous variants get +0.15;
   - every pair among the top 3 heterozygous variants is emitted on **one row** with score mean + 0.15;
   - a lone het gets −0.10, so it ranks below any recessive-compatible row.
7. **EPCR.** Row scores are scaled to (0, 1] (top = 0.95) with a small monotone decrement, so every row has a distinct threshold. A maximum of 10 rows is allowed; here only 3 candidate rows existed.

## Methods-form answers

- **Generative AI:** Anthropic Claude Code (model: Claude Opus 5.5) was used interactively to review the plan, write and debug the pipeline code, and draft this report. Plan/tier: Claude Max subscription (consumer terms), used through the Claude Code CLI. Data handling: the session ran on the local machine; the assistant did not upload any dataset files. The assistant ran the pipeline locally. Panel-region variant annotations were therefore visible to the assistant during the session, but no raw reads and no whole-genome VCF content were sent to it. No patient data was sent to any external annotation service: all annotation ran locally.
- **Automated or curated:** The submission file is the **unedited automated output** of `mva_track1.py rank`. The only manual inputs are the gene list and the scoring weights, both fixed before the final run.
- **Manual review:** After ranking, the two *BUB1B* variants were checked by hand for read support (AD, GQ, FILTER), for gnomAD and ClinVar consistency, and for agreement with the phenotype. Nothing was changed as a result.
- **Data:** Public data only.
  - Ensembl release 113 (GRCh38 FASTA and GTF; REST API for gene coordinates)
  - ClinVar GRCh38 VCF (2026-10-19)
  - gnomAD v4.1 exome and genome sites
  - AlphaMissense hg38 (Cheng et al., *Science* 2023)
  - OMIM gene–disease entries for MVA1–3
  - Ensembl VEP 113 (Docker image `ensemblorg/ensembl-vep:release_113.0`)
  - bcftools/htslib/samtools 1.24
- **Proprietary data:** none.
- **Compound heterozygotes:** Yes. The method emits proposed compound-het pairs on a single row (`chrom_2..alt_2`). Homozygous variants and single heterozygous variants are emitted as single-variant rows.
- **Secondary or incidental findings:** None were reported. The analysis is restricted to the MVA gene panel, so incidental findings outside it were not searched for. All rows are `primary`.
- **Run time and cost:** Building the reference bundle (`scripts/prepare_reference.sh`) takes about 10 minutes, almost all of it download time. The pipeline itself takes under 2 minutes on an 8 GB Apple-silicon laptop (VEP under Docker/x86 emulation). It needs about 250 MB of reference data and has no cloud compute cost.

## Abstract (≤500 words)

We ranked candidate causal variants for a proband with clinically suspected mosaic variegated aneuploidy (MVA). The pipeline is a transparent, reproducible gene-panel approach driven by a recessive-inheritance model. The single-sample GATK VCF (GRCh38) is restricted to the three established MVA genes (*BUB1B*, *CEP57*, *TRIP13*) and four related mitotic-checkpoint and kinetochore genes, with ±5 kb padding. Variants are normalised with bcftools and annotated locally with Ensembl VEP 113, using ClinVar, gnomAD v4.1 exome and genome frequencies, and AlphaMissense. To avoid the 25 GB VEP cache and to keep patient data on the local machine, every annotation resource is sliced to the panel regions by remote tabix queries. This makes the whole resource bundle about 250 MB and the full run a few minutes on a laptop.

Variants pass genotype QC (PASS, DP ≥ 10, GQ ≥ 20), a population-frequency ceiling of 1% (maximum over gnomAD genetic ancestry groups), removal of ClinVar benign variants, and a consequence filter. Surviving variants receive an interpretable weighted score combining predicted severity (LoF, or missense scaled by AlphaMissense), rarity, ClinVar evidence and a gene prior. Because MVA is autosomal recessive, variants are then grouped by gene. Pairs of heterozygous variants in the same gene are reported together as one compound-heterozygous row, and lone heterozygous variants are down-weighted. Scores are mapped to strictly decreasing EPCR values.

Of 252 panel variants, two passed all filters, and both are in *BUB1B*. One is a ClinVar pathogenic/likely-pathogenic nonsense variant, c.2210T>G p.(Leu737Ter), very rare in gnomAD. The other is a novel or ultra-rare kinase-domain missense variant, c.3006T>G p.(Asn1002Lys), with AlphaMissense 0.92. Both are high-quality heterozygous calls with balanced allele fractions. Their combination ranks first. The recorded phenotype (an embryonal-type soft-tissue cancer, growth restriction, prematurity) is consistent with *BUB1B*-related MVA1. The phenotype was used only to confirm the genotype-derived ranking, not to drive it.

**Strengths:** the pipeline is fully automated and reproducible from one script, and every score component is interpretable. Compound-het pairs are modelled explicitly, as the scoring requires. It uses public data only, no patient data leaves the machine, and it runs on a laptop with minimal compute.

**Limitations:**
- **Phase is unknown.** The two variants are about 11 kb apart and the single-sample VCF has no read-backed phasing for them. Trans configuration is assumed, not shown; parental testing or long-read phasing would be needed.
- **Narrow panel.** Genes outside the panel are not assessed.
- **No splicing or non-coding analysis.** SpliceAI was not available, so deep-intronic, UTR and promoter variants (relevant for hypomorphic *BUB1B* alleles) were filtered out unless they were coding.
- **No structural variants or CNVs.** These were not assessed.
- **Phenotype not used in scoring.** Exomiser-style phenotype matching was not used because of local disk limits.
- **Hand-set weights.** The score weights are expert-set, not learned.
