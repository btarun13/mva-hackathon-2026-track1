# Plan – Pathway enrichment, protein consequences and STRING network for the MVA proband
**Rare Disease, Real Kid: MVA Hackathon 2026** · drafted 2026-10-23 · builds on `track1_plan.md` and `work/vep.tsv`

## Goal
Explain, at the pathway and protein level, how the proband's variants change the cell compared with a typical person: which proteins are altered, which complexes and processes they sit in, which interactions are expected to be lost or weakened, and what the predicted functional outcome is (here, loss of spindle-assembly-checkpoint signalling). Then ask whether any of the affected proteins or their partners are druggable.

## Mentor guidance folded into this plan (2026-10-23)
- Convert variant gene IDs to official gene names with HGNC (https://www.genenames.org/), then explore **STRING**, **GSEA**, **DGIdb** (https://dgidb.org/) and **ChEMBL**.
- **Variants that change a protein are the most important input.** The key analysis is a comparison of the *expected* protein-interaction network when no variant changes a protein (a typical person) against the network with this child's altered proteins. Everything below is organised around that comparison.

## Starting point and a caveat that shapes the whole plan
`work/vep.tsv` covers only the 7-gene candidate panel ±5 kb. It yields 11 unique genes (CENPE, KIF18A, BUB1B, CENATAC, BUB1, CEP57, TRIP13, plus neighbours BDH2, TRAPPC4, FAM76B and one unnamed non-coding gene) and only 2 variants that change a protein in a meaningful way, both in *BUB1B*:

| Variant | Protein effect | Expected consequence |
|---|---|---|
| chr15:40209701 T>G, c.2210T>G | p.(Leu737Ter), premature stop in exon 17/23 | Transcript likely degraded by nonsense-mediated decay; any escaping protein lacks residues 737–1050, i.e. the whole C-terminal kinase domain |
| chr15:40220612 T>G, c.3006T>G | p.(Asn1002Lys), exon 23/23 | Missense inside the kinase domain, AlphaMissense 0.92; predicted to impair kinase activity or fold |

**Enrichment on a hand-picked panel is circular.** Any gene-set test on those 11 genes will "discover" mitotic checkpoint and kinetochore terms because that is how the panel was chosen. The plan therefore has two tiers:

- **Tier A (genome-wide, meaningful statistics):** re-run annotation on the *whole* VCF, derive an unbiased list of genes carrying rare, protein-altering variants, and test that list for enrichment. This answers "is the mitotic-checkpoint signal there without us putting it in?" and finds any second hit or modifier.
- **Tier B (panel, mechanistic):** take the two *BUB1B* variants and ask what breaks downstream, using protein domains, STRING interaction partners, and pathway membership. This answers "what is happening in this child's cells?"

Both tiers compare the proband to a "normal person" by using **gnomAD as the population baseline** (frequency and constraint) rather than a matched control sample, because no control genome is provided.

## Tier A – Genome-wide rare protein-altering gene set

1. **Annotate the whole VCF.** Reuse `mva_track1.py annotate` but without the panel restriction. This needs the full VEP cache (~25 GB) or a cloud VEP run; the cache-free GTF mode will not scale. Keep `--pick --canonical`, ClinVar, gnomAD v4.1, AlphaMissense; add LOFTEE if available.
2. **Filter to a "burden" gene list.** Per variant: PASS, DP ≥ 10, GQ ≥ 20; MAX_AF < 0.1 % (stricter than Track 1 because the genome-wide list must stay short); consequence HIGH or MODERATE; drop ClinVar benign. Per gene, flag genes that are **recessive-compatible** (homozygous or ≥ 2 heterozygous qualifying variants) separately from genes with a single heterozygous hit. Expect roughly 100–300 genes with any hit and 5–30 recessive-compatible genes.
3. **Background.** Use all protein-coding genes covered by the VCF at adequate depth, not the whole genome, to avoid coverage bias.
4. **Over-representation analysis (ORA).** Run g:Profiler (`gprofiler-official`) or `gseapy.enrichr` offline with GO Biological Process, Reactome, KEGG, and the HPO gene–phenotype set. Test three lists: all hit genes, recessive-compatible genes, and the Track 1 panel as a positive control. Benjamini–Hochberg FDR < 0.05.
5. **Rank-based GSEA (optional).** Build a per-gene score (max variant score from `variant_score`, 0 for untouched genes) and run pre-ranked GSEA with `gseapy.prerank`. This uses the whole genome as background and avoids an arbitrary cutoff.
6. **Expected and informative outcomes.** If *BUB1B* is the only recessive-compatible gene in checkpoint-related terms, that supports a single-gene diagnosis. If other recessive-compatible genes fall into mitosis, DNA repair, or ciliary terms, record them as possible modifiers or a differential diagnosis.

## Tier B – Protein and network consequences of the *BUB1B* variants

0. **Gene-name mapping (HGNC).** Map every Ensembl gene ID in `work/vep.tsv` (and later the genome-wide table) to the approved HGNC symbol, previous symbols and UniProt accession via the genenames.org REST API (`/fetch/ensembl_gene_id/<ENSG>`). Keep a `results/hgnc_map.tsv`. This makes all downstream STRING, enrichment, DGIdb and ChEMBL queries use one consistent identifier and resolves unnamed genes such as ENSG00000289202.

1. **Protein mapping.** Use UniProt O60566 (BUBR1, 1050 aa) and InterPro/Pfam domains: N-terminal TPR/KEN boxes (Cdc20 binding, MCC formation), KARD (PP2A-B56 binding), C-terminal kinase domain (~766–1050). Mark both variants on this map. Check AlphaFold DB (AF-O60566) for the local structure around Asn1002, and gnomAD constraint (pLI, LOEUF, missense o/e) to state how intolerant *BUB1B* is to loss of function.
2. **Predicted protein states.** Allele 1: no protein (NMD) or a truncated 736-residue protein that still has the MCC-forming N-terminus but no kinase domain and no KARD. Allele 2: full-length protein with a kinase-domain substitution. Note that BUBR1's kinase domain is considered a pseudokinase in humans, so the pathogenic mechanism of Asn1002Lys is more likely destabilisation and reduced protein level, which fits the hypomorphic mechanism reported for MVA1 missense alleles. Cite the literature for MVA1 compound-het genotypes (one null + one hypomorphic allele is the classic pattern).
3. **STRING network.** Query the STRING API (v12, *Homo sapiens*, score ≥ 0.7, experimental+database evidence) with:
   - the seed set: BUB1B, BUB1, BUB3, MAD1L1, MAD2L1, CDC20, TRIP13, CEP57, CENPE, KIF18A, CENATAC, TTK, PLK1, AURKB, KNL1, PPP2R5B/C/E;
   - the Tier A recessive-compatible genes, to see whether any touch this module.
   Retrieve interaction partners and STRING's own functional enrichment for the network. Export as TSV and a PNG via the API's `network` image endpoint; also save the Cytoscape-compatible edge list.
4. **Complex-level interpretation.** Annotate which STRING partners are in the mitotic checkpoint complex (MCC: BUBR1–BUB3–MAD2–CDC20), the KNL1–BUB1–BUBR1 kinetochore recruitment axis, the PP2A-B56 chromosome-alignment axis, and the TRIP13–p31comet MCC-disassembly axis. State which of these the two variants should disrupt: loss of BUBR1 abolishes MCC formation and the anaphase-wait signal, and weakens kinetochore–microtubule error correction via PP2A-B56.
5. **Normal vs proband interactome (the central comparison).** Build two networks from the same STRING query:
   - **Reference network ("normal person"):** all seed proteins present and intact, all STRING edges kept. This is the expected interactome when no variant changes a protein.
   - **Proband network:** mark proteins carrying protein-changing variants (from `work/vep.tsv`: BUB1B only; from Tier A: any genome-wide rare protein-altering hits). For each altered protein, classify its edges by the interacting region using UniProt/IntAct binding-site annotations and the literature: edges mediated by regions that are lost (for the Leu737Ter allele: the KARD and the kinase domain, so PP2A-B56 and kinase-domain partners) are scored *lost*; edges mediated by retained regions (TPR/KEN boxes: CDC20, BUB3, MAD2, KNL1) are scored *weakened* because protein level is reduced by NMD and the missense allele is hypomorphic.
   - **Diff table:** for every edge touching an altered protein, report partner, STRING score, evidence channels, the mediating region, and the predicted state (intact / weakened / lost), plus the downstream complex affected. Also list first-neighbour proteins whose complexes become incomplete (MCC, PP2A-B56 at kinetochores).
   - Figures: side-by-side network PNGs with lost edges dashed and altered nodes coloured; exported edge lists for Cytoscape.
6. **Druggability (DGIdb + ChEMBL).** Query the DGIdb GraphQL API for BUB1B and every first-neighbour partner in the diff table (CDC20, BUB1, BUB3, MAD2L1, TTK, PLK1, AURKB, PPP2CA/PPP2R5x, TRIP13, CENPE, KIF18A). Pull ChEMBL target records and the number of approved/clinical compounds for each. Interpret with care: MVA is a loss-of-function disease, so most hits will be checkpoint *inhibitors* used in oncology (e.g. TTK, PLK1, AURKB, CENPE, KIF18A inhibitors) which would worsen the defect in normal cells but are relevant to treating the aneuploid tumour; flag anything that could *stabilise* residual BUBR1 or compensate for it. Output `results/druggability.tsv`.
7. **Normal vs proband narrative.** Produce a two-column summary: in a typical cell, unattached kinetochores recruit BUBR1 into the MCC, inhibit APC/C-CDC20, and anaphase waits; in the proband, BUBR1 is absent or hypomorphic, the MCC is weak, anaphase proceeds with unattached chromosomes, and daughter cells gain or lose chromosomes (the variegated aneuploidy). Link this to the phenotype items: aneuploidy-driven embryonal cancer, growth restriction, and the parents' recurrent miscarriages (carrier status does not explain miscarriage, so flag that as unexplained).

## Deliverables
- `results/hgnc_map.tsv` (Ensembl ID → HGNC symbol → UniProt)
- `work/genome_vep.tsv` and `work/burden_genes.tsv` (Tier A; patient-derived, not committed)
- `results/enrichment_{all,recessive,panel}.tsv` and dot plots (committed: gene lists in aggregate only)
- `results/string_edges_reference.tsv`, `results/string_edges_proband.tsv`, `results/interactome_diff.tsv`, side-by-side network PNGs, `results/string_enrichment.tsv`
- `results/BUB1B_protein_map.png` (domain diagram with the two variants and the partner-binding regions)
- `results/druggability.tsv` (DGIdb + ChEMBL for altered proteins and first neighbours)
- `pathway_string_report.md`: methods, figures, the normal-vs-proband table, limitations

## Tooling
- Python via `uv`: `pandas`, `gseapy`, `gprofiler-official`, `requests` (HGNC, STRING, DGIdb GraphQL, ChEMBL REST APIs), `networkx`, `matplotlib`
- VEP 113 with the GRCh38 cache (Docker), or a cloud VEP run for Tier A
- Optional: Cytoscape with the stringApp for a publication-grade network figure

## Limitations to state up front
- No control genome; gnomAD is the baseline, so population-stratification effects are possible.
- Tier B is hypothesis-level: no RNA or protein data to confirm NMD or reduced BUBR1 levels.
- Phase of the two *BUB1B* variants is still unknown; trans configuration is assumed.
- STRING edges are population-level knowledge, not patient-specific measurements; the network shows what *could* be affected, not what *is*.
- Enrichment of the panel genes (Tier B positive control) is by construction and must not be reported as a finding.

## Order of work
1. Tier B steps 1–2 (a day, no new data needed) → gives the mechanistic core even if Tier A stalls.
2. STRING query and figure (half a day).
3. Tier A genome-wide annotation (depends on cache download or cloud time) → enrichment → fold into the report.
4. DGIdb/ChEMBL druggability pass over the diff-table proteins.
5. Write `pathway_string_report.md`.
