# MVA Hackathon 2026 — Track 1 report (all-variants companion file)

**Team / user:** Tarun2712
**Model 1:** `vep-rules-v1` (full-table variant: submission file `Tarun2712_vep-rules-v1_all-variants.csv`)
**Code:** https://github.com/btarun13/mva-hackathon-2026-track1 (`scripts/run_pipeline.sh` builds `work/vep.tsv`; `scripts/all_vep_to_submission.py` then writes this file byte-for-byte)
**Relation to the main submission:** `submissions/Tarun2712_vep-rules-v1.csv` is the filtered, 3-row file intended for scoring. This file is the same pipeline **without any filtering or row cap**: every one of the 252 panel variants annotated by VEP is emitted as one row, in the same 12-column format, with the same scoring weights. It exists so that the complete ranked table behind the main submission can be inspected.

## Result

| Rank | EPCR | Variant (GRCh38) | Gene | Zygosity | Interpretation |
|---|---|---|---|---|---|
| 1 | 0.950 | chr15:40209701 T>G | *BUB1B* | het | c.2210T>G p.(Leu737Ter), stop_gained, ClinVar P/LP, gnomAD max AF 1.0e-4 |
| 2 | 0.764 | chr15:40220612 T>G | *BUB1B* | het | c.3006T>G p.(Asn1002Lys), missense, AlphaMissense 0.92, gnomAD max AF 9.0e-7 |
| 3 | 0.492 | chr15:40216470 A>G | *BUB1B* | het | c.2679-1026A>G, deep intronic, absent from gnomAD; dropped by the main pipeline (no SpliceAI support) |
| 4 | 0.394 | chr15:40192892 C>T | *BUB1B* | het | c.1059-3653C>T, deep intronic, ClinVar VUS, max AF 7.5e-3; dropped by the main pipeline |
| 5 | 0.350 | chr2:110671177 G>A | *BUB1* | het | c.423-609C>T, intronic, max AF 8.8e-5; dropped by the main pipeline |
| 6 | 0.310 | chr2:110668376 C>T | *BUB1* | het | c.568-527G>A, intronic, max AF 2.1e-4; dropped by the main pipeline |
| 7 | 0.309 | chr4:103119480 G>T | *CENPE* | het | c.7329+668C>A, intronic, max AF 3.3e-4; dropped by the main pipeline |
| 8–252 | 0.294 → 0.001 | … | all panel genes | het/hom | common intronic, UTR, upstream/downstream and synonymous variants; 34 ClinVar benign |

Only ranks 1 and 2 survive the main pipeline's filters. They are the two *BUB1B* variants that form the compound-heterozygous candidate in the main submission. In this file they appear as single-variant rows, because every variant is listed individually and no pair rows are generated.

**Composition of the 252 rows**

| Breakdown | Count |
|---|---|
| Heterozygous / homozygous-alt | 173 / 79 |
| Consequence class: LoF / missense / splice-region / other (intronic, UTR, flanking, synonymous) | 1 / 3 / 2 / 246 |
| By gene: *CENPE* / *KIF18A* / *BUB1B* / *CENATAC* / *BUB1* / *CEP57* / *TRIP13* / neighbouring genes (*BDH2*, *TRAPPC4*, *FAM76B*) / no gene symbol | 136 / 48 / 17 / 16 / 14 / 4 / 3 / 13 / 1 |
| ClinVar: P/LP / VUS / benign or likely benign / not reported | 1 / 1 / 35 / 215 |
| FILTER ≠ PASS (QD2, MQ40, LowQual) | 3 |
| Would be removed by the main pipeline | 250 |

Reasons the main pipeline removes a variant, as recorded in the `notes` column (a variant can have several):

| Reason | Rows |
|---|---|
| Non-coding or synonymous without SpliceAI support | 247 |
| gnomAD max AF > 1 % | 242 |
| ClinVar benign / likely benign | 35 |
| VCF FILTER not PASS | 3 |

## Approach

Steps 1–3 (gene panel, VCF preparation, VEP 113 annotation) are identical to the main submission and are not repeated here; see `submissions/Tarun2712_track1_report.md`. This file differs from step 4 onward.

4. **No filters.** QC, rarity, ClinVar and consequence filters are evaluated but **not applied**. Instead, each failed filter is written into the `notes` column after the text `pipeline would drop:` so that the reader can see exactly why the main submission excludes the variant.
5. **Variant score.** Same formula as the main submission: score = 0.35·severity + 0.20·rarity + 0.20·ClinVar + 0.25·gene prior. Two differences are needed so that no variant is dropped:
   - ClinVar benign variants receive a ClinVar contribution of 0.0 instead of being removed.
   - Non-coding and synonymous variants keep their severity of 0.6·SpliceAI (0.0 here, since SpliceAI was unavailable) instead of being removed.
   Rarity for variants above 1 % is 0.0, as before. The gene prior is 1.0 for *BUB1B*, *CEP57*, *TRIP13*; 0.5 for *CENPE*, *KIF18A*, *CENATAC*, *BUB1*; 0.1 for neighbouring genes caught by the ±5 kb padding.
6. **No recessive grouping.** Every variant is a single-variant row (`chrom_2..alt_2` empty). The homozygous bonus, compound-het pair rows and lone-het penalty of the main pipeline are not applied, since the aim here is a complete per-variant table rather than a diagnosis.
7. **EPCR.** Scores are sorted descending (ties broken by chromosome and position), scaled so the top row is 0.95, then decremented by 0.001 per rank and floored at 0.001. Every row therefore has a distinct, strictly decreasing EPCR, and all 252 values lie in (0, 1].

## Methods-form answers

- **Generative AI:** Anthropic Claude Code (Claude Opus 5.5 for the pipeline and the main report; Claude Fable 5.1 for this companion script and report) was used interactively, under a Claude Max subscription through the Claude Code CLI, with the default consumer data setting ("Help improve Claude" enabled). All processing ran on the local machine. Panel-region annotations and the resulting table were visible to the assistant; no raw reads or whole-genome VCF content were sent to it, and no patient data was sent to any external annotation service.
- **Automated or curated:** The file is the **unedited automated output** of `scripts/all_vep_to_submission.py`. No rows were added, removed or reordered by hand.
- **Manual review:** The two *BUB1B* coding variants at ranks 1–2 were reviewed by hand for the main submission. The remaining 250 rows were not individually reviewed; the `notes` column records the automated reason for their low ranking.
- **Data:** Public data only, identical to the main submission (Ensembl 113, ClinVar 2026-10-19, gnomAD v4.1 exomes and genomes, AlphaMissense hg38, OMIM MVA1–3, Ensembl VEP 113 Docker image, bcftools/htslib/samtools 1.24).
- **Proprietary data:** none.
- **Compound heterozygotes:** Not emitted in this file. The compound-het pair is represented in the main submission; here its two component variants occupy ranks 1 and 2 as separate rows.
- **Secondary or incidental findings:** None reported. All rows are `primary`. Analysis is restricted to the MVA gene panel plus genes overlapping the ±5 kb padding.
- **Run time and cost:** After `scripts/run_pipeline.sh` (under 2 minutes), the companion script runs in a few seconds on a laptop. No cloud compute cost.
- **Validity for the official scorer:** The file follows the required column layout but has 252 rows, which exceeds the 10-row limit enforced by `mva_track1.py validate` and the official scorer. It is a companion table for inspection, **not** a file intended for upload.

## Abstract (≤500 words)

This companion file reports the complete per-variant table behind our Track 1 submission for a proband with clinically suspected mosaic variegated aneuploidy (MVA). The main submission applies a recessive-aware, rules-based filter to a candidate-gene panel (*BUB1B*, *CEP57*, *TRIP13*, plus *CENPE*, *KIF18A*, *CENATAC*, *BUB1*) and reports three rows. Here the same pipeline is run with every filter switched off, so that all 252 panel variants annotated by Ensembl VEP 113 appear as individual rows in the submission format, ranked by the same interpretable score (0.35 severity, 0.20 rarity, 0.20 ClinVar, 0.25 gene prior) and mapped to strictly decreasing EPCR values.

The table confirms that the main result is not an artefact of aggressive filtering. The two highest-scoring variants are the *BUB1B* nonsense variant c.2210T>G p.(Leu737Ter), ClinVar pathogenic/likely pathogenic and very rare in gnomAD, and the *BUB1B* missense variant c.3006T>G p.(Asn1002Lys), ultra-rare with AlphaMissense 0.92. They outrank every other panel variant by a wide margin (EPCR 0.950 and 0.764 against 0.492 for rank 3). The next five rows are rare deep-intronic variants in *BUB1B*, *BUB1* and *CENPE* that the main pipeline excludes because no SpliceAI evidence was available. Below them, 242 of the 252 variants exceed 1 % population frequency in at least one gnomAD ancestry group, 35 are ClinVar benign or likely benign, and 3 fail the VCF quality filter. The distribution of genes (136 *CENPE*, 48 *KIF18A*, 17 *BUB1B*) simply reflects gene length and the padding used for region extraction.

**Strengths:** complete transparency. Every variant in the panel is visible with its score and, where applicable, the exact reason the main pipeline removes it. The file is produced by a short, deterministic script from the same intermediate VEP table as the main submission, so the two files are guaranteed consistent.

**Limitations:**
- **Not a scorable submission.** The row count exceeds the 10-row limit; the file is for review only.
- **No compound-het rows.** Pair rows are deliberately omitted, so the recessive model of the main submission is not expressed here.
- **Non-coding variants are under-scored.** Without SpliceAI, every intronic, UTR and synonymous variant receives zero severity, so the two rare deep-intronic *BUB1B* variants at ranks 3–4 cannot be distinguished from genuinely neutral variants.
- **Same narrow panel, no structural variants, hand-set weights** as the main submission.
