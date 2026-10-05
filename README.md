# hedgehog-variability

Single-cell RNA-seq analysis of why genetically identical cells respond
differently to the same signal.

This repository contains the analysis notebooks for:

> Ramu A, Cohen BA. *Transcription factor fluctuations underlie cell-to-cell
> variability in a signaling pathway response.* GENETICS (2023).
> Preprint: [bioRxiv doi:10.1101/2022.11.30.518555](https://www.biorxiv.org/content/10.1101/2022.11.30.518555v1)

<!-- TODO(Avi): add the GENETICS journal DOI link once confirmed. -->

## The question

Clonal cells grown in the same environment still behave differently. Some
respond quickly to a signal and others slowly or not at all. We asked whether
random fluctuations in **transcription factor (TF) activity** explain part of
this variability.

As a model we used the **Hedgehog (HH) signaling pathway** in NIH3T3-CG mouse
fibroblasts, which carry a GFP reporter of pathway activity. Cells were treated
with the Smoothened agonist SAG and profiled by flow cytometry and 10x Genomics
single-cell RNA-seq over a time course.

## Key findings

- NIH3T3-CG cells contain **fast-responding and slow-responding substates**
  with distinct expression profiles.
- Overexpressing the candidate TFs Prrx1, Snai1, Jun or Srebf2 each misexpressed
  a significant number of fast-responder genes, and overexpressing **Prrx1**
  alone was enough to **speed up the HH response**.
- Fluctuations in TF activity are a plausible general source of phenotypic
  variability among genetically identical cells.

## Analysis approach

| Step | Methods |
| --- | --- |
| QC, normalization, cell-cycle scoring, clustering, UMAP | [Seurat](https://satijalab.org/seurat/) |
| Differential expression between responder states | Seurat, hypergeometric enrichment tests |
| Trajectory / pseudotime analysis | [Monocle 3](https://cole-trapnell-lab.github.io/monocle3/), [Slingshot](https://bioconductor.org/packages/slingshot/), [TSCAN](https://bioconductor.org/packages/TSCAN/) |
| Supervised cell-state classification | [Garnett](https://cole-trapnell-lab.github.io/garnett/), trained on labeled fast and slow responders and applied to new experiments |

## Notebooks

| File | What it does |
| --- | --- |
| [`analyze.Rmd`](analyze.Rmd) | Analysis of the first scRNA-seq experiment (SAG for 0, 17 and 30 hours): QC, clustering, cell-cycle filtering and trajectory analysis |
| [`mark_competent_cells.Rmd`](mark_competent_cells.Rmd) | Labels cells as fast or slow responders and finds differentially expressed genes between them, with pathway enrichment |
| [`run_garnett_classifier.Rmd`](run_garnett_classifier.Rmd) | Trains a Garnett classifier on fast and slow responders and uses it to predict the state of cells in the Prrx1 experiment |
| [`analyze_prrx1_combined.Rmd`](analyze_prrx1_combined.Rmd) | Analysis of cells with and without doxycycline-induced Prrx1 expression at 0, 19 and 30 hours of SAG |
| [`run_slingshot.R`](run_slingshot.R) | Pseudotime ordering of cells with Slingshot |

## Requirements

R with these packages: Seurat, monocle3, garnett, slingshot, TSCAN, scran,
SingleCellExperiment, org.Mm.eg.db, biomaRt, dplyr, ggplot2, RColorBrewer.

## Data

The notebooks were run on a computing cluster and read intermediate results
from absolute paths (for example `/scratch/bclab/...`). To re-run them, update
those paths to point to your local copy of the data.

<!-- TODO(Avi): add the GEO accession for the raw and processed scRNA-seq data. -->
