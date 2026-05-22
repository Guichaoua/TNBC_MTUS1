# Transcriptomic profile of MTUS1-low TNBC reveals candidate therapeutic strategies

This repository contains the analysis code, selected notebooks, derived data tables, and figures for the article:

**Transcriptomic profile of MTUS1-low TNBC reveals candidate therapeutic strategies**

Authors: Gwenn Guichaoua, Olivier Collier, Sylvie Rodrigues-Ferreira, Clara Nahmias, and Véronique Stoven.

The study is framed as computational hypothesis generation. Tumour transcriptomic analyses define a reproducible MTUS1-low programme enriched for MYC-linked and stress-adaptation pathways. DepMap CRISPR dependency data are then used as an external functional prioritisation layer. Previously published WEE1/PKMYT1 biology is interpreted as external biological support, not as the primary discovery of this analysis.

## Repository layout

- `scripts/`: publication-facing scripts for figure regeneration.
- `src/`: reusable analysis helpers from the original project.
- `notebooks/`: selected exploratory notebooks retained for auditability.
- `data/derived/`: small derived tables used in figures and supplementary tables.
- `data/raw/`: placeholder for downloaded public raw data; raw data are not committed. DepMap raw files should be placed under `data/raw/Cell-lines/`.
- `figures/`: final manuscript figures and selected source figure inputs.

## Quick reproduction

Create the environment:

```bash
conda env create -f environment.yml
conda activate mtus1-tnbc-article
```

Regenerate article figures from included notebooks and raw DepMap files:

```bash
python scripts/regenerate_article_figures.py
python scripts/recompose_article_figure1.py
```

The `regenerate_article_figures.py` script and the reusable DepMap helpers expect raw DepMap files under `data/raw/Cell-lines/`.
Required files for the dependency figures include `CRISPRGeneEffect.csv`,
`OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected.csv`, and
`res/meta_24cellLines_chronos.csv`.

Download and preprocess TCGA cohorts from the NCI Genomic Data Commons:

```bash
python src/TCGA_load.py TCGA-BRCA
python src/TCGA_load.py TCGA-BRCA,TCGA-OV,TCGA-PRAD
```

This creates `data/<TCGA-project>/treated_files/counts.txt`,
`data/<TCGA-project>/treated_files/genes.txt`, and
`data/<TCGA-project>/treated_files/metadata.txt`. Intermediate GDC downloads
are kept under `data/<TCGA-project>/download_data/` and are not intended to be
committed. The same workflow is documented in `notebooks/TCGA_load.ipynb`.

## Data availability

All raw datasets used in the manuscript are public. GEO datasets are available under GSE181466, GSE192341, and GSE202203. TCGA data are available through the NCI Genomic Data Commons. SRA-derived datasets SRP042620 and SRP157974 were accessed through recount3. DepMap transcriptomic and Chronos CRISPR-Cas9 dependency data were obtained from the DepMap portal.

The repository includes derived tables sufficient to audit the main plotted results and supplementary dependency tables. Large raw expression matrices and downloaded TCGA/GEO/SRA files are intentionally excluded.

## Citation notes

When submitting, cite the exact DepMap release used in the manuscript and include the DepMap portal URL. The current manuscript cites DepMap Public 24Q4; if analyses are refreshed to another release, update the release citation in the manuscript bibliography and Methods.

## Manuscript status

This cleaned repository is intended for peer-review transparency and does not include the article source. The Overleaf-ready manuscript package is maintained separately. The original analyses were developed in notebooks, with publication figures consolidated into scripts for reproducibility.
