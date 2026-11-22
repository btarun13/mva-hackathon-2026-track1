# MVA proband: pathway, protein-interaction and druggability analysis
**Rare Disease, Real Kid: MVA Hackathon 2026 — follow-up to Track 1** · Tarun2712 · 2026-10-24
**Plan:** `pathway_string_plan.md` · **Code:** `scripts/tier_a_genomewide.sh`, `scripts/gnomad_fetch_chunks.sh`, `scripts/add_gnomad_maxaf.py`, `scripts/tier_a_analysis.py`, `scripts/tier_b_interactome.py` · **Outputs:** `results/`

## Summary

The question was: at the level of proteins, complexes and pathways, what is different in this child's cells compared with a typical person, and is any of it druggable?

- **Genome-wide (Tier A):** of 104,972 coding-region variants, 176 are rare (gnomAD max ancestry-group AF < 0.1 %) and change a protein. They fall in 164 genes; only 10 genes are compatible with recessive inheritance (homozygous/hemizygous, or ≥ 2 heterozygous hits). *BUB1B* is the top-ranked of the 10 and the only one with a loss-of-function allele plus a ClinVar pathogenic entry. The other nine are explained by known artefact classes (HLA/KIR polymorphism, adjacent variants on one haplotype, a low-quality call) or lack any phenotype link. No pathway is significantly enriched among the 164 genes after multiple-testing correction, which is the expected result for a single-gene recessive disease. The genome-wide scan therefore confirms the Track 1 answer without having been told where to look.
- **Protein level (Tier B):** the proband carries exactly two rare protein-changing variants in the mitotic-checkpoint module, both in BUBR1 (*BUB1B*): p.Leu737Ter, which removes the C-terminal pseudokinase domain and is expected to trigger nonsense-mediated decay, and p.Asn1002Lys, inside that domain (AlphaMissense 0.92). All 21 high-confidence STRING interactions of BUBR1 are predicted to be weakened by reduced protein dosage; the BUBR1–CENPE interaction, which depends on the pseudokinase domain, is predicted lost. The complexes affected are the mitotic checkpoint complex (BUBR1–BUB3–MAD2–CDC20), its inhibition of APC/C^CDC20, the KNL1–BUB1–BUBR1 kinetochore axis, PP2A-B56 recruitment via the KARD, and CENPE-dependent chromosome congression.
- **Normal vs proband:** in a typical cell an unattached kinetochore recruits BUBR1 into the MCC, APC/C^CDC20 is held off and anaphase waits. In the proband BUBR1 is absent from one allele and destabilised on the other, so the wait signal is weak, anaphase proceeds with mis-attached chromosomes, and daughter cells gain or lose whole chromosomes: mosaic variegated aneuploidy.
- **Druggability:** BUBR1 itself has no potent compounds in ChEMBL and one weak DGIdb entry. Its partners TTK/MPS1, PLK1, AURKB, CDK1, CENPE and KIF18A have clinical-stage *inhibitors*. These oppose the checkpoint and would aggravate the germline defect in normal cells, but are the relevant class for the aneuploid tumour; nothing in DGIdb or ChEMBL restores BUBR1 function.

## 1. Inputs and methods

**Patient data.** Single-sample GRCh38 GATK VCF (5,012,204 records). Restricted to protein-coding exons ± 25 bp from Ensembl 113 (223,561 merged intervals, 95 Mb) → 104,015 records; contigs renamed to `chr*`; multiallelics split and indels left-aligned with bcftools → 104,972 variant records. Phenotype used only for interpretation, never for ranking.

**Annotation (cache-free VEP 113).** The 25 GB VEP cache could not be downloaded at the available bandwidth, so VEP ran in GTF mode with the full Ensembl 113 GTF and primary-assembly FASTA (`--pick --canonical --hgvs`), ClinVar (GRCh38 VCF, 2026-10-23) as a custom track and AlphaMissense as a plugin. gnomAD v4.1 exome and genome frequencies (AF, AF_grpmax) were fetched for the 11,979 HIGH/MODERATE positions only, by remote tabix range queries against the public Google-storage files, in 100-region chunks with retries (an earlier single-stream fetch silently truncated on 27 of 48 chromosome files; the chunked fetch completed 161/161 chunks). MAX_AF = max(AF, AF_grpmax) over both sets; 12,003 of 12,015 HIGH/MODERATE variants were found in gnomAD.

**Filters and scoring (Tier A).** FILTER PASS, DP ≥ 10, GQ ≥ 20, het or hom-alt, protein-coding transcript; IMPACT HIGH or MODERATE; MAX_AF < 0.001 (absent from gnomAD counts as rarest); ClinVar benign/likely-benign removed. Each variant gets the Track 1 score without the gene prior (0.35 severity + 0.20 rarity + 0.20 ClinVar). Genes are *recessive-compatible* if they carry ≥ 1 homozygous/hemizygous or ≥ 2 heterozygous qualifying variants.

**Enrichment.** Over-representation (hypergeometric, Benjamini–Hochberg) with gseapy against GO Biological Process 2023, Reactome 2022, KEGG 2021 and the Human Phenotype Ontology, using as background the 15,529 genes that carry at least one QC-pass coding variant in this genome (coverage-aware background). Pre-ranked GSEA was attempted and dropped: 97 % of background genes tie at score 0, so the ranking is arbitrary and the result uninterpretable. The 7-gene Track 1 panel was also run through Enrichr as a positive control; because the panel was chosen for checkpoint function, that result is circular and is reported only as a sanity check.

**Identifiers.** Every Ensembl gene ID was mapped to its HGNC symbol and UniProt accession through the genenames.org REST API (`results/hgnc_map.tsv`); all genome-wide hits already carried a symbol.

**Interactome (Tier B).** STRING v12.5 functional network, score ≥ 0.7, for a seed module of 24 checkpoint/kinetochore proteins. BUBR1 regions and the partners they mediate were taken from UniProt O60566 and the literature (KEN1, TPR, D-box, KEN2, GLEBS, ABBA, KARD, pseudokinase). Each STRING edge of an altered protein was classified as *lost* if its mediating region is removed by the truncation or hit by the missense, otherwise *weakened* by dosage. Druggability from the DGIdb GraphQL API and the ChEMBL REST API (targets, potent activities with pChEMBL ≥ 6, mechanism compounds and max clinical phase).

## 2. Tier A — genome-wide rare protein-altering variants

### 2.1 Funnel

| Stage | Variants | Genes |
|---|---|---|
| Coding-region records after normalisation | 104,972 | |
| QC-pass, protein-coding transcript | 99,877 | 15,529 (background) |
| HIGH or MODERATE impact | 11,029 | |
| gnomAD max AF < 0.1 %, not ClinVar benign | **176** | **164** |
| Recessive-compatible genes | | **10** |

### 2.2 The ten recessive-compatible genes

| Rank | Gene | Hits (zyg.) | Best variant(s) | gnomAD max AF | Assessment |
|---|---|---|---|---|---|
| 1 | ***BUB1B*** | 2 het | p.Leu737Ter (ClinVar P/LP) + p.Asn1002Lys (AlphaMissense 0.92) | 1.0e-4 / 9.0e-7 | **Causal candidate.** Only gene with LoF + ClinVar pathogenic; matches MVA1 |
| 2 | *ADAMTS1* | 2 het | two frameshifts at codons 220 and 225, identical AF | 8.8e-5 | Adjacent indels with identical frequency: one complex event split into two records, likely restoring frame. No MVA link |
| 3 | *SERPINA1* | 7 het | cluster of 7 variants within codons 312–347, AFs from AMR/AFR genomes only | ≤ 5.6e-4 | Dense cluster on one haplotype in a region with known paralogous/alignment issues; artefactual |
| 4 | *MRPL23* | 1 hom | p.Arg162Ser, absent from gnomAD | absent | Mitochondrial ribosomal protein; no AlphaMissense score, no ClinVar; no phenotype link |
| 5 | *DNAH10* | 2 het | two missense, AF 5e-4 and 8.5e-4 | – | Very large axonemal dynein; near the frequency cut-off; no phenotype link |
| 6 | *HLA-DQA1* | 2 het | two missense | 6.5e-4 / 2.3e-4 | Hyper-polymorphic locus; mapping artefact class |
| 7 | *ZNF619* | 2 het | p.Gln556Lys + p.Asp557Val, adjacent, identical AF | 9.3e-4 | One multi-nucleotide variant on one haplotype |
| 8 | *IGBP1* | 1 hemi (X) | p.Met151Leu, ClinVar VUS, AlphaMissense 0.16 | 1.0e-4 | X-linked; proband is male (X het fraction 0.15, Y calls present). IGBP1/α4 binds PP2A-C and is the only non-panel gene with a STRING edge into the checkpoint module (IGBP1–PPP2CA 0.999). Predicted benign; noted as a possible modifier only |
| 9 | *ZNF418* | 1 hom | p.Ser622Arg, DP 10, GQ 33 | 9.2e-4 | Low-quality call at the QC floor; disregard |
| 10 | *KIR3DL3* | 2 het | two missense | 4.3e-4 / 9.0e-7 | Hyper-polymorphic KIR locus; mapping artefact class |

Full table: `results/tierA_gene_burden.tsv` (gene level). Per-variant table with genotypes: `work/genome/rare_protein_altering_variants.tsv` (patient-level, not committed).

### 2.3 Enrichment

- **All 164 genes with a rare protein-altering variant:** no term passes FDR < 0.05 in any library. The closest are HPO skin-fragility terms driven by *ITGB4/LAMA3/COL7A1* (FDR 0.058) and the KEGG immune terms driven by HLA/KIR (FDR 0.06). The only checkpoint-related term, GO "centromeric sister chromatid cohesion" (*BUB1B*, *SGO1*), has FDR 0.12. This is the expected picture for a monogenic recessive disorder: one gene, no pathway-wide burden.
- **Ten recessive-compatible genes:** 71 nominally significant terms, but they decompose into (i) KEGG immune terms carried by *HLA-DQA1* and *KIR3DL3*, the artefact pair, and (ii) HPO terms each carried by *BUB1B* alone (nephroblastoma, rhabdomyosarcoma, increased nuchal translucency, leiomyosarcoma). Single-gene overlaps are not pathway evidence; they simply restate that *BUB1B* is an MVA gene.
- **Panel positive control (circular):** the 7 panel genes return "spindle assembly checkpoint signalling" (GO, FDR 1.3e-6), "mitotic spindle checkpoint" (Reactome, FDR 2.8e-7) and HPO "nephroblastoma" (FDR 7.8e-7), confirming the libraries and the pipeline behave as expected. Not a finding.

Figure: `results/tierA_enrichment_ora.png`. Tables: `results/tierA_enrichment_ora.tsv`, `results/enrichment_panel_positive_control.tsv`.

### 2.4 STRING overlay and druggability of the recessive-compatible genes
Only *IGBP1* connects to the checkpoint seed module (via PPP2CA). DGIdb lists drugs for 6 of the 10 genes, almost all generic (HLA-DQA1 pharmacogenomic associations, alpha-1-antitrypsin for *SERPINA1*); none is relevant to the mechanism. `results/tierA_string_links_to_checkpoint_module.tsv`, `results/tierA_dgidb_recessive_genes.tsv`.

## 3. Tier B — what the two *BUB1B* variants do to the protein and its network

### 3.1 Protein-changing variants in the panel
Of 252 panel variants, four change a protein. Two are common benign *CENPE* polymorphisms (p.Thr2090Met, p.Phe1535Leu; AF 0.20, ClinVar benign). The other two are the *BUB1B* alleles below; they are the only altered proteins in the module. `results/protein_changing_variants.tsv`.

| Allele | Change | Where on BUBR1 (1050 aa) | Predicted protein state |
|---|---|---|---|
| 1 | c.2210T>G p.(Leu737Ter), exon 17/23 | between the KARD (665–682) and the pseudokinase domain (766–1050) | Premature stop > 50 nt upstream of the last exon junction → nonsense-mediated decay expected; any escaping product lacks the pseudokinase domain |
| 2 | c.3006T>G p.(Asn1002Lys), exon 23/23 | inside the pseudokinase domain | Full-length but destabilised; human BUBR1's kinase domain is catalytically inactive and its main role is to stabilise the protein, so missense here lowers BUBR1 levels (the classic hypomorphic MVA1 allele) |

Figure: `results/BUB1B_protein_map.png`.

### 3.2 Reference vs proband interactome

STRING seed module: 23 proteins, 128 edges at score ≥ 0.7. Its own enrichment (sanity check) returns "mitotic spindle checkpoint" (Reactome, FDR 2.5e-33), "kinetochore" (GO CC, FDR 1.9e-18), KEGG "cell cycle" (FDR 6.6e-34) and the DISEASES term "mosaic variegated aneuploidy syndrome" (3/3 genes, FDR 3.9e-6). `results/string_enrichment.tsv`.

Of BUBR1's 21 edges:

| Predicted state in proband | n | Partners | Mediating region |
|---|---|---|---|
| **Lost / severely impaired** | 1 | CENPE | pseudokinase domain (removed by allele 1, mutated in allele 2) |
| Weakened by dosage, binding region retained | 14 | CDC20, MAD2L1, BUB3, KNL1, BUB1, ANAPC1, CDC27, PLK1, CDK1, TTK, PPP2CA, PPP2R5B/C/E | KEN/ABBA/D-box, TPR, GLEBS, KARD |
| Weakened by dosage, region not annotated | 6 | MAD1L1, AURKB, TRIP13, MAD2L1BP, KIF18A, ZW10 | – |

Figure: `results/interactome_reference_vs_proband.png` (left: reference network; right: proband, lost edge dotted red, weakened edges dashed orange). Edge table: `results/interactome_diff.tsv`; reference edges: `results/string_edges_reference.tsv`.

### 3.3 What is different from a normal person

| Step in mitosis | Typical cell (reference interactome) | Proband (BUBR1 null + hypomorph) |
|---|---|---|
| Kinetochore recruitment | KNL1 MELT motifs → BUB3–BUB1 → BUBR1 via its TPR/GLEBS regions | Binding regions intact but BUBR1 supply is reduced; fewer BUBR1 molecules at unattached kinetochores |
| Checkpoint signal | BUBR1 + BUB3 + MAD2 + CDC20 assemble the MCC; BUBR1's KEN/ABBA/D-box motifs block APC/C^CDC20 | MCC forms less efficiently; APC/C^CDC20 inhibition is leaky; the "wait anaphase" signal is weak |
| Error correction | KARD-bound PP2A-B56 counteracts Aurora B, stabilising correct attachments | Less KARD-bound PP2A at kinetochores → persistent attachment errors |
| Congression | BUBR1 pseudokinase domain binds and localises CENPE | Interaction lost → impaired chromosome alignment |
| Checkpoint silencing | TRIP13/p31^comet disassemble the MCC once attached | Unchanged, but acts on a smaller MCC pool, so silencing is faster |
| Outcome | Anaphase only after all kinetochores attach; euploid daughters | Anaphase with unattached/mis-attached chromosomes → random gains and losses → mosaic variegated aneuploidy |

Clinical fit: aneuploidy-driven embryonal tumour (rhabdomyosarcoma-type, HP:0002859), pre- and postnatal growth restriction and prematurity are the *BUB1B*-MVA1 pattern. Parental carrier status does not by itself explain the recorded recurrent miscarriages; that remains unexplained.

### 3.4 Druggability (DGIdb + ChEMBL)

| Protein | Role relative to BUBR1 | ChEMBL potent activities / max phase | DGIdb drugs (approved) | Relevance |
|---|---|---|---|---|
| **BUB1B** | altered | 0 / – | 1 (0) | No tool compound or drug; not directly druggable |
| AURKB | KARD/PP2A axis antagonist | 2,412 / 3 | 99 (11) | Inhibitors in oncology trials |
| TTK (MPS1) | phosphorylates KARD | 2,066 / 2 | 17 (0) | Checkpoint inhibitors (BAY 1217389, BOS172722) |
| PLK1 | phosphorylates KARD | 1,118 / 3 | 191 (19) | Volasertib etc. |
| CDK1 | phosphorylates KARD | 1,201 / 3 | 97 (7) | Pan-CDK inhibitors |
| KIF18A | congression motor | 1,066 / – | 0 | Clinical inhibitors for aneuploid tumours (not yet in DGIdb) |
| CENPE | lost partner | 80 / – | 2 (0) | GSK-923295 |
| PPP2CA | KARD-bound phosphatase | 14 / 2 | 7 (2) | LB-100 (PP2A inhibitor) |
| CDC20, MAD2L1, BUB3, KNL1, TRIP13, MAD1L1 | MCC / kinetochore | ≤ 4 / – | 0–2 | Undruggable today |

Interpretation: every actionable target in this module is a *checkpoint kinase inhibitor*. In the child's normal tissues that direction is wrong (it would further weaken an already weak checkpoint). In the tumour, which is likely highly aneuploid and checkpoint-dependent, inhibitors of KIF18A, TTK, PLK1 or AURKB are the mechanistically coherent class, and KIF18A inhibitors are being developed specifically for chromosomally unstable cancers. No agent stabilises BUBR1 or compensates for its loss. `results/druggability.tsv`, `results/dgidb_interactions.tsv`.

## 4. Limitations
- **Coding regions only.** Tier A looked at exons ± 25 bp; deep-intronic, UTR, promoter, structural and copy-number variation were not assessed genome-wide.
- **Frequencies for HIGH/MODERATE variants only**; LOW/MODIFIER variants have no gnomAD annotation in this run (not needed for the analysis).
- **Phase unknown.** The two *BUB1B* alleles are assumed in trans. Parental testing or read-backed phasing is required to prove it.
- **No functional data.** NMD, reduced BUBR1 level and loss of CENPE binding are predictions from sequence, domain annotation and the literature, not measurements in this patient.
- **STRING edges are population knowledge.** The diff table shows which interactions *should* be affected given where the variants fall; it is not a patient-specific proteomics result.
- **Enrichment power.** With one causal gene, over-representation tests cannot be significant and are reported as negative controls; the recessive-list "hits" are artefact-driven.
- **gnomAD baseline, no control genome.** "Normal" is defined by population frequency and annotated reference biology.

## 5. Reproducibility
```bash
scripts/run_pipeline.sh                       # Track 1 (panel) -> work/vep.tsv
uv run --no-project --with pandas,requests,networkx,matplotlib python scripts/tier_b_interactome.py
scripts/tier_a_genomewide.sh                  # coding-region VEP (GTF mode) + chunked gnomAD lookup
uv run --no-project --with pandas,requests,gseapy,networkx,matplotlib python scripts/tier_a_analysis.py
```
Large references live on the USB drive `/Volumes/SAN_DRIVE/mva_refs` (FASTA, GTF, ClinVar, AlphaMissense; ~3 GB). Patient-derived per-variant tables stay in `work/` and are git-ignored. Generative AI: Claude Code (Claude Fable 5.1) was used interactively for planning, coding and drafting; all computation ran locally except the public STRING, HGNC, DGIdb and ChEMBL API queries, which received gene names only.
