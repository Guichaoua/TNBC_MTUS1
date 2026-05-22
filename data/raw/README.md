Raw public datasets are not committed to this repository because the full analysis workspace is large.

Place downloaded raw files here when re-running the complete analysis. The figure-generation script expects DepMap files under:

`data/raw/Cell-lines/`

Required DepMap files for the dependency figures include:

- `CRISPRGeneEffect.csv`
- `OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected.csv`
- `res/meta_24cellLines_chronos.csv`

Reusable helper functions in `src/Depmap_utils.py` also expect optional
DepMap/PRISM/GDSC files in this same `data/raw/Cell-lines/` folder when those
legacy analyses are rerun.

TCGA cohorts are downloaded through the GDC API with `src/TCGA_load.py`. The
processed files used by the cohort notebooks are written under:

`data/<TCGA-project>/treated_files/`

For example:

```bash
python src/TCGA_load.py TCGA-BRCA
```

This creates `counts.txt`, `genes.txt`, and `metadata.txt` for the requested
project. Intermediate downloaded GDC files are placed under
`data/<TCGA-project>/download_data/`.

The manuscript cites DepMap Public 24Q4. If a newer release is used, update both the manuscript methods and bibliography.
