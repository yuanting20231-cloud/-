# Microbiome figure workflows

This repository contains reproducible R code for two microbiome figure workflows based on the uploaded metadata, ASV abundance matrix, and taxonomy table.

## Available figures

### Figure 2

1. Chao1 alpha-diversity boxplots for rhizosphere and root.
2. Bray–Curtis PCoA with 95% confidence ellipses and PERMANOVA.
3. Family-level relative-abundance stacked bars.
4. Rhizosphere family-level differential abundance.
5. Root family-level differential abundance.

Run with:

```r
source("install_packages.R")
source("run_all.R")
```

### Figure 3

The new workflow reproduces the three analytical modules in the supplied network figure:

1. **Figure 3a** — network connectivity after random removal of increasing proportions of nodes, shown separately for rhizosphere and root.
2. **Figure 3b** — four NetShift-style circular comparison networks for `R vs S` and `R vs E` in rhizosphere and root. Edge colours distinguish control-only, case-only, and shared associations. Node size represents the Neighbor Shift score; highlighted nodes have increased scaled betweenness and high NESH values.
3. **Figure 3c** — UpSet-style intersections of driver families among the four NetShift comparisons.

Run with:

```r
source("install_packages.R")
source("run_figure3.R")
```

The workflow writes separate PNG/PDF panels, a combined figure, network tables, model statistics, driver-family tables, and `sessionInfo()` to `output/`.

## Uploaded files checked

- Metadata: 18 samples in six groups (`RSV`, `RRV`, `SSV`, `SRV`, `ESV`, and `ERV`), with three replicates per group.
- Abundance matrix: 24,122 ASVs × 18 samples, no missing or negative values, and 31,556 reads per sample.
- The abundance matrix is stored losslessly in four compact archive parts and is reconstructed automatically on first use.

## Required input files

```text
data/
├── sample_metadata.tsv
├── taxonomy.tsv
├── asv_counts.compact.b64.part001
├── asv_counts.compact.b64.part002
├── asv_counts.compact.b64.part003
└── asv_counts.compact.b64.part004
```

### Metadata interpretation

The current working interpretation is:

- second character `S` = `Rhizosphere`
- second character `R` = `Root`
- first character `R`, `S`, or `E` = treatment identity

Edit `data/sample_metadata.tsv` if these biological meanings differ.

### Taxonomy

Place the uploaded taxonomy table at:

```text
data/taxonomy.tsv
```

The first column must be `#ID`, and the table must contain a `family` column. Figure 3 is calculated at family level.

### ASV abundance matrix

`R/00_restore_abundance.R` reconstructs the following file automatically:

```text
data/asv_abundance.tsv
```

No manual decompression or format conversion is required.

## Figure 3 code structure

```text
R/
├── 00_config.R
├── 00_helpers.R
├── 06_network_helpers.R
├── 07_fig3a_network_robustness.R
├── 08_fig3b_netshift.R
└── 09_fig3c_driver_upset.R
run_figure3.R
```

## Figure 3 default comparisons

`R` is currently treated as the focal/case group and is compared with `S` and `E` as reference/control groups within each compartment. Edit `NETSHIFT_COMPARISONS` in `R/00_config.R` to reverse or rename these groups.

## Figure 3 methods

- Counts are aggregated to family level.
- Families are filtered by prevalence and mean relative abundance.
- Counts are CLR transformed after adding a pseudocount.
- Because each treatment-compartment group currently contains only three samples, networks are estimated with shrinkage correlations rather than ordinary correlations.
- The strongest associations are retained using an adaptive absolute-correlation threshold and edge-count limits defined in `R/00_config.R`.
- Network robustness is measured as the number of retained edges relative to the maximum possible edges in the original network after random node removal.
- The NetShift-style NESH score uses neighborhood overlap and directional enrichment of case-specific neighbours. Driver families require a positive change in scaled betweenness and `NESH >= 1.0` by default.
- The UpSet panel summarizes intersections among the four driver-family sets.

## Important statistical limitation

The current design has only three replicates per treatment-compartment group. The code will run using regularized correlations, but group-specific network topology, NESH scores, driver-family calls, confidence bands, and regression statistics should be treated as exploratory. For publication-level network inference, substantially more independent samples per group are recommended. The thresholds in `R/00_config.R` should then be re-evaluated, and a compositional network method with resampling or stability selection should be used as a sensitivity analysis.
