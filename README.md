# Microbiome Figure 2: five-panel R workflow

This repository contains reproducible R code for the five graph modules in the supplied reference figure:

1. **Panel a** — Chao1 alpha-diversity boxplots, separated into rhizosphere and root.
2. **Panel b** — Bray–Curtis PCoA with 95% confidence ellipses and PERMANOVA statistics.
3. **Panel c** — Family-level relative-abundance stacked bars using the top 20 families plus `Others`.
4. **Panel d** — Rhizosphere family-level differential abundance, with group means, mean differences, 95% confidence intervals, and Student's t-test P values.
5. **Panel e** — Root family-level differential abundance using the same statistical layout.

## Uploaded files checked

- The sample grouping file contains 18 samples and six groups: `RSV`, `RRV`, `SSV`, `SRV`, `ESV`, and `ERV`, with three replicates per group.
- The two previously uploaded large text files are identical taxonomy tables containing `#ID`, phylum, class, order, family, genus, and species columns.
- The newly uploaded ASV abundance matrix has been added to this repository. It contains **24,122 ASVs × 18 samples**, has no missing or negative values, and each sample contains **31,556 reads**.

## Repository structure

```text
microbiome_fig2_code/
├── data/
│   ├── sample_metadata.tsv
│   ├── taxonomy.tsv                         # copy the taxonomy table here
│   ├── asv_counts.compact.b64.part001       # uploaded abundance matrix archive
│   ├── asv_counts.compact.b64.part002
│   ├── asv_counts.compact.b64.part003
│   └── asv_counts.compact.b64.part004
├── R/
│   ├── 00_config.R
│   ├── 00_restore_abundance.R               # automatically creates asv_abundance.tsv
│   ├── 00_helpers.R
│   ├── 01_panel_a_chao1.R
│   ├── 02_panel_b_pcoa.R
│   ├── 03_panel_c_family_composition.R
│   ├── 04_panel_d_rhizosphere_differential.R
│   └── 05_panel_e_root.R
├── install_packages.R
└── run_all.R
```

## Input data

### 1. Metadata

`data/sample_metadata.tsv` has already been prepared from the uploaded grouping file.

The current working interpretation is:

- second character `S` = `Rhizosphere`
- second character `R` = `Root`
- first character `R`, `S`, or `E` = treatment identity

This is an explicit working assumption because the grouping file contains only `Sample` and `Group`. Correct `Treatment` or `Compartment` directly in `sample_metadata.tsv` when the biological meanings differ.

### 2. Taxonomy

Copy either one of the two identical uploaded taxonomy files to:

```text
data/taxonomy.tsv
```

The first column must be `#ID`, and the file must contain a `family` column. Panels c–e require this file.

### 3. ASV abundance matrix

The uploaded abundance matrix is stored losslessly in four compact archive parts. When `R/00_config.R` is sourced, `R/00_restore_abundance.R` automatically reconstructs:

```text
data/asv_abundance.tsv
```

The restored table has ASVs in rows, samples in columns, and non-negative integer read counts as cell values. No manual conversion is required.

## Run

From the repository root:

```r
source("install_packages.R")
source("run_all.R")
```

The first run reconstructs `data/asv_abundance.tsv` automatically. All PNG and PDF outputs are written to `output/`.

## Comparisons for panels d and e

The default comparisons are:

- `R - S`
- `R - E`

Edit `COMPARISONS` in `R/00_config.R` if another treatment should be the reference. The plotted difference is always `group1 - group2`.

## Statistical methods

- Panel a: Chao1 estimated with `vegan::estimateR`; Dunn post hoc grouping letters with BH correction.
- Panel b: Bray–Curtis dissimilarity, Lingoes-corrected PCoA, and PERMANOVA with 999 permutations.
- Panels d/e: family-level relative abundance; equal-variance Student's t-test; unadjusted P values displayed to reproduce the reference layout; BH-adjusted values retained in the internal results table.

## Important note

With only three replicates per treatment, 95% ellipses and family-level t tests can be unstable. Interpret P values and confidence intervals cautiously, and consider compositional differential-abundance methods such as ANCOM-BC2 or ALDEx2 as a sensitivity analysis for publication.
