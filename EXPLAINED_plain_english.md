# What we did, explained from scratch

This is a plain-language walkthrough of the whole Track 1 analysis: what the problem was, what each step did, what we found, and (importantly) what we did *not* do, what we guessed at, and what is still uncertain. Nothing here is softened. If something is a weakness, it is called a weakness.

---

## 1. The problem in one paragraph

A child has a rare disease called **mosaic variegated aneuploidy (MVA)**. In MVA, the machinery that checks chromosomes are shared out correctly when a cell divides is faulty, so many of the child's cells end up with the wrong number of chromosomes. The hackathon gave us the child's **genome sequence** and a short description of their **symptoms**, and asked one question: **which specific DNA change (or pair of changes) is causing the disease?** We had to submit a ranked list of guesses, and the organisers score it against the answer the clinicians already confirmed.

---

## 2. What we were given

| Input | What it is | Did we use it? |
|---|---|---|
| `WGS_EX2312012_HGWCNDSX7.vcf.gz` (315 MB) | A **VCF file**: a list of every place where the child's DNA differs from the standard human reference genome. About **5 million** differences. Each line says: chromosome, position, what the reference has, what the child has, and whether the child has it on one copy or both copies of that chromosome. | **Yes, this was the main input.** |
| 8 FASTQ files (~85 GB) | The raw sequencing reads, i.e. the data the VCF was made from. | **No.** Not needed; we got a clear answer from the VCF. |
| `Challenge_Clinical_Phenotype_1.docx` | The symptoms: a soft-tissue cancer in childhood, kidney calcium deposits, short stature, failure to thrive, premature birth, very low birth weight, and recurrent miscarriages in the parents. | **Only to sanity-check the answer at the end.** It was not fed into the scoring. |
| The disease name: MVA | | **Yes, this drove the whole strategy.** |

Everyone carries millions of DNA differences and almost all of them are harmless. The job is to find the one or two that matter. 5 million candidates is far too many to look at, so the real work is narrowing down.

---

## 3. The strategy: don't search the whole genome, search where the answer almost certainly is

Two facts about MVA shaped everything:

1. **Only a handful of genes are known to cause it.** Medical databases (OMIM) list three: **BUB1B**, **CEP57** and **TRIP13**. There are a few more genes with weaker evidence or related biology (we added CENPE, KIF18A, CENATAC, BUB1).
2. **MVA is recessive.** You need *both* copies of the gene broken to get the disease. That can happen two ways:
   - the same change on both copies (**homozygous**), typical when parents are related; or
   - two *different* damaging changes, one on each copy (**compound heterozygous**). This is the common situation.

So instead of 5 million candidates, we looked only inside those 7 genes. That shrinks the problem to **252 DNA changes**. And because the disease is recessive, we specifically looked for genes with *two* damaging changes, not one.

**Why this is a reasonable bet:** the disease was already diagnosed clinically, so the answer is very likely in a known MVA gene. **Why it's a limitation:** if the real cause were in some gene nobody has linked to MVA yet, our approach could never find it. We accepted that trade-off.

---

## 4. The steps, one by one

All the logic lives in one Python script, `mva_track1.py`. `scripts/run_pipeline.sh` runs the steps in order.

### Step 0: Get the tools and the reference data

- **bcftools / samtools / tabix**: standard command-line tools for cutting, filtering and reading genome files. Installed with Homebrew.
- **VEP (Variant Effect Predictor)**: a program from Ensembl that, for each DNA change, works out *what it does* to a gene (e.g. "this changes amino acid number 737 into a STOP signal"). Normally VEP needs a 25 GB database. Your laptop only had 13 GB free, so instead we ran VEP inside Docker and gave it only the pieces it needed: the gene map (GTF) and DNA sequence (FASTA) for the five chromosomes that carry our seven genes.
- **Annotation databases**, all public, each cut down to just our gene regions (about 250 MB total instead of many gigabytes):
  - **gnomAD v4.1**: how common each DNA change is in ~800,000 people. A change seen in 5% of healthy people cannot cause a rare recessive disease on its own.
  - **ClinVar**: a database where clinical labs record whether a given change has already been seen in patients and classified as "pathogenic", "benign", etc.
  - **AlphaMissense**: a DeepMind AI model that predicts, for amino-acid changes, how likely they are to be damaging (score 0 to 1).

### Step 1: Gene list → genome coordinates

`genes.txt` lists the seven genes. The script asks Ensembl's public web service where each gene sits on the genome and writes a small file of regions, padded by 5,000 bases on each side so that control regions next to the gene are included too.

(Only the gene *names* were sent to Ensembl. Nothing about the patient was.)

### Step 2: Cut the patient's VCF down to those regions

- The patient's file names chromosomes as `15`; the scoring system and our reference files use `chr15`. We renamed them.
- `bcftools norm` tidied up the formatting: one record that listed two alternative alleles on one line was split into two lines, and one small deletion was shifted to its leftmost possible position (standard practice so the same change is always written the same way).
- Result: **252 DNA changes** inside the seven genes, with the child's genotype for each (one copy or two copies, how many reads supported it, and the quality score).

### Step 3: Annotate every change with VEP

For each of the 252 changes VEP added: which gene, which transcript, what the effect is (intronic, missense, stop-gained, etc.), its gnomAD frequency, its ClinVar status, and its AlphaMissense score. Result: a big table, `work/vep.tsv`.

### Step 4 & 5: Filter

The script then threw out changes that are very unlikely to be causal:

| Filter | Rule | Why | Left after |
|---|---|---|---|
| (start) | | | 252 |
| Quality | the variant caller must have flagged it PASS | 3 were flagged low quality | 249 |
| Read support | at least 10 reads covering it, genotype quality ≥ 20 | avoid sequencing noise | 244 |
| Rarity | gnomAD frequency under 1% in every population | a recessive-disease allele this severe cannot be common | 10 |
| Effect | must actually alter the protein (stop, frameshift, missense, splice site) | intronic/non-coding changes usually do nothing | **2** |
| ClinVar | drop anything ClinVar calls benign | (nothing was dropped by this) | 2 |

The rarity filter did almost all the work: 234 of the 244 good-quality changes are common in the population, so they're just normal human variation.

### Step 6: Score, pair up, rank

The two survivors were given a score (0 to 1) combining four things: how severe the predicted effect is (35%), how rare it is (20%), ClinVar evidence (20%), and how strongly the gene is linked to MVA (25%).

Then the recessive logic: both survivors are in **the same gene (BUB1B)** and the child has **one copy of each**. That is exactly the compound-heterozygous pattern we were looking for. The script put the two together on one row as a pair (with a bonus for being recessive-compatible), and also listed each one alone as a fallback.

Finally, scores were converted to the "EPCR" number the organisers want (just a confidence between 0 and 1, highest first).

---

## 5. What we found

**Two changes in the gene BUB1B, one on each copy (most likely):**

| | Change 1 | Change 2 |
|---|---|---|
| Position (GRCh38) | chr15:40,209,701 T→G | chr15:40,220,612 T→G |
| In the protein | **p.Leu737Ter**: amino acid 737 becomes a STOP signal, so the protein is cut off early and the whole kinase (enzyme) end is missing. | **p.Asn1002Lys**: amino acid 1002 changes from asparagine to lysine, inside the kinase domain. |
| ClinVar | Already classified **Pathogenic / Likely pathogenic** by clinical labs. | Never reported. |
| How common | ~1 in 13,000 people carry it (gnomAD). | ~1 in 1.5 million (seen once in gnomAD exomes). Essentially novel. |
| AlphaMissense | (n/a for a stop) | **0.92** = "likely pathogenic" |
| Read support | 46 reads, 21 ref / 25 alt, quality 99 | 28 reads, 15 ref / 13 alt, quality 99 |

**Why this is convincing:**
- BUB1B is *the* classic MVA gene (MVA type 1).
- The combination fits the known pattern for BUB1B-MVA: one allele that completely kills the protein, plus a second, milder allele.
- The specific cancer in the symptom document is the one most characteristic of BUB1B-related MVA, as opposed to the other two MVA genes.
- Both calls are clean: roughly half the reads carry the change, which is what you expect for a change on one of two copies.
- Nothing else survived the filters. It isn't "the best of a messy list"; it's the only candidate.

**What we can't prove from this data:** that the two changes are on *different* copies of chromosome 15 (called being "in trans"). If both were on the same copy, the other copy would be fine and this would not explain the disease. They're 11,000 bases apart and the sequencing reads are too short to link them, so the VCF can't tell us. Testing the parents (each should carry one) would settle it. We've flagged this clearly in the report. It is the single biggest scientific caveat.

---

## 6. What's in the submission

| File | What it is |
|---|---|
| `submissions/Tarun2712_vep-rules-v1.csv` | The prediction. Row 1: the BUB1B pair, confidence 0.95. Rows 2–3: each variant alone, as backup. This is the raw, unedited script output. |
| `submissions/Tarun2712_track1_report.md` | The write-up judges will read. |
| `submissions/Tarun2712_methods_description_form.xlsx` | The organisers' question form, filled in with the same answers. |
| A GitHub repo (not created yet) | Code so others can reproduce the result. A local git repo with one commit is ready to push. |

**How scoring works:** rank 1 correct = 100 points. A compound-het pair only gets full credit if **both** variants are on the **same row**; that's why row 1 is a pair. Rows 2 and 3 are there so that if the organisers' answer is only one of the two variants, we still get half credit.

---

## 7. The honest list

### Things I changed in the plan/code to make it run
1. **VEP without its database.** The plan assumed the 25 GB VEP cache. Your disk couldn't hold it, so I added a "cache-free" mode that feeds VEP a gene map + DNA sequence + sliced databases. You approved this. Side effect: the "maximum population frequency" number is computed slightly differently from VEP's built-in one (gnomAD's `AF_grpmax` ignores a few very small populations). For these two variants it makes no difference, both are extremely rare either way.
2. **One bug fix in `subset`.** The original command would always crash (`bcftools view -R` can't read from a pipe). Changed `-R` to `-T`. You approved this.
3. **Chromosome renaming step.** Not in the plan; needed because the patient file says `15` and everything else says `chr15`.

### Things in the plan that were NOT done
- **SpliceAI** (predicts whether a change breaks RNA splicing). The precomputed score files require an Illumina account. Skipped.
- **Exomiser** (ranks genes by how well the symptoms match). Its database wouldn't fit on the disk. Skipped. The symptoms were only used by me, by eye, to check that BUB1B makes sense.
- **Plan B** (re-process the 85 GB of raw reads). Not needed; the VCF gave a clean answer.

### Known weaknesses I pointed out at the start and deliberately left alone (you asked me to)
- **The filters would discard a non-coding second allele.** In some BUB1B families the second allele is a change in a regulatory region, not in the protein-coding part. Our filter removes those. It turned out not to matter here because a coding second allele was found, but if it hadn't been, the pipeline would have missed it.
- **`--pick` plus 5 kb padding** can occasionally assign a variant to a neighbouring gene (you can see BDH2, FAM76B and TRAPPC4 in the output; these are neighbours, not panel genes). Didn't affect the result.
- **The SpliceAI parser in the script reads the wrong columns.** Irrelevant this run because SpliceAI wasn't used, but it would be wrong if someone turns it on.
- **Hard quality cut-offs** (DP ≥ 10, GQ ≥ 20) on such a small panel could in principle drop a true variant. Here the 8 variants removed by quality were all common/non-coding anyway.
- **The score weights (35/20/20/25) are hand-picked**, not learned from data. They're reasonable but arbitrary. With only two candidates, the weights didn't change the outcome.

### Things that are uncertain or that I could not verify
- **Whether the answer is right.** It is very well supported, but only the organisers' scorer can confirm it. I've not submitted anything.
- **Phasing** (see section 5). We assume the two variants are on different copies. The data cannot show this.
- **Whether the scorer's truth uses exactly our coordinates.** The template uses `chr15` and GRCh38, and our rows match that format. The one indel we left-aligned isn't in the answer, so normalisation can't have changed the coordinates that matter.

### Things you still need to do
1. **Fill in the AI-use line** in the report (search for `[TODO`). The organisers require it. It must state the Claude plan/tier you are on and the data-handling setting. I don't know these, so I didn't invent them.
2. **Create the GitHub repo and push**, then put the URL in the report header. The local repo is committed and `data/`/`work/` are git-ignored so the patient's data can't be pushed by accident. The repo has to go public after the hackathon ends.
3. **Revoke the Hugging Face token you pasted into the chat** and make a new one. It's in the conversation transcript.
4. **Delete the patient data** when you're done, and email the organisers, as you agreed when you accepted the dataset terms.
5. Read the report yourself before submitting. It's written as your work. You should be comfortable defending every sentence.

### About AI involvement, plainly
Claude (an AI assistant, me) reviewed the plan, downloaded the public reference data, wrote the scripts and helper code, ran the pipeline, inspected the results, and wrote the report, the form and this document. You made the decisions at each fork (which VEP route, whether to fix the bug, using the token). The patient's variant annotations for the 252 panel positions were visible to me during the session; the raw reads and the full genome file were not sent anywhere. The organisers explicitly allow AI assistance as long as it's disclosed, which the report does.

---

## 8. Small glossary

- **Genome / reference genome**: a person's full DNA; the "reference" is a standard agreed sequence everyone is compared against. GRCh38 is the current version.
- **Variant**: any place where someone's DNA differs from the reference.
- **VCF**: the file format listing variants.
- **Heterozygous (het)**: the change is on one of your two copies of a chromosome. **Homozygous (hom)**: on both.
- **Compound heterozygous**: two *different* changes in the same gene, one on each copy. Together they knock the gene out.
- **Recessive**: a disease that needs both copies of a gene to be broken.
- **Missense**: a change that swaps one amino acid for another. May or may not matter.
- **Stop-gained / nonsense**: a change that inserts a premature STOP; the protein is truncated. Usually severe.
- **Intronic / non-coding**: changes outside the parts of the gene that encode protein. Usually harmless.
- **gnomAD**: a huge database of variants seen in the general population. Used to tell "normal variation" from "rare".
- **ClinVar**: a database of clinicians' verdicts on specific variants.
- **AlphaMissense**: an AI predictor of whether a missense change is damaging.
- **VEP**: the tool that figures out what each variant does to a gene.
- **Phase / in trans / in cis**: whether two variants are on different copies (trans, needed for compound het) or the same copy (cis).
- **EPCR**: the organisers' name for "your estimated probability that this is the cause", 0 to 1.
