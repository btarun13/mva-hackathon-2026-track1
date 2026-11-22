# MVA Hackathon 2026 — Track 1 (variant prediction)

Candidate-gene, recessive-aware variant ranking for the MVA proband. Single script: `mva_track1.py`
(subcommands `genes`, `subset`, `annotate`, `rank`, `validate`). Method details: `submissions/Tarun2712_track1_report.md`.

## Requirements
- bcftools / htslib / samtools ≥ 1.20 (`brew install bcftools htslib samtools`)
- Docker (VEP 113: `docker pull --platform linux/amd64 ensemblorg/ensembl-vep:release_113.0`; wrapper in `bin/vep`)
- [uv](https://docs.astral.sh/uv/) (pulls pandas on the fly)
- Access to the gated dataset `SageBio/mva-hackathon-2026-data` (HF login)

## Reproduce
```bash
hf download SageBio/mva-hackathon-2026-data WGS_EX2312012_HGWCNDSX7.vcf.gz WGS_EX2312012_HGWCNDSX7.vcf.gz.tbi \
  --repo-type dataset --local-dir ./data
scripts/prepare_reference.sh          # ~10 min, builds ref/ (~250 MB, public data only)
scripts/run_pipeline.sh <hf-username> # <2 min, writes submissions/<hf-username>_vep-rules-v1.csv
```

## Layout
| Path | What |
|---|---|
| `mva_track1.py` | pipeline (normalise/subset, VEP annotation, filter, score, compound-het grouping, EPCR, validation) |
| `genes.txt` | candidate gene panel (core MVA genes + related checkpoint genes) |
| `scripts/prepare_reference.sh` | panel-sliced Ensembl/ClinVar/gnomAD/AlphaMissense bundle (no 25 GB VEP cache needed) |
| `scripts/run_pipeline.sh` | end-to-end run |
| `bin/vep` | Docker wrapper for VEP 113 |
| `track1_plan.md` | original analysis plan |

Patient data (`data/`, `work/`) is git-ignored and must not be committed.
