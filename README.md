# CAVD-CD74-MMP9

Analysis code for the manuscript:

**Integrated bioinformatics and machine learning identify CD74 and MMP9 as candidate immune-related tissue biomarkers of calcific aortic valve disease**

Cong Huo, Hong Zuo, Xiaohui Shi, Wei Zhou, Wenbing Hou, Gang Fan (corresponding author)

Third Department of Cardiology, Xianyang Central Hospital, Xianyang, Shaanxi, China

## Overview

The pipeline integrates bulk transcriptomics, machine learning, immune deconvolution, single-cell RNA-seq and exploratory Mendelian randomisation (MR) to nominate CD74 and MMP9 as candidate tissue-level immune-related biomarkers of calcific aortic valve disease (CAVD).

## Data sources (all public)

| Dataset | Platform | Use |
|---|---|---|
| GSE83453 | Microarray (GPL10558), 19 calcified vs 8 normal aortic valves | Discovery |
| GSE153555 | RNA-seq, aortic valve (AS / AI / NC) | Cross-platform replication & model transferability |
| GSE12644 | Microarray (GPL570), 10 calcified vs 10 normal | Fully untouched third cohort |
| GSE220774 | Single-cell RNA-seq, 82,356 cells, 5 donor valves | Cellular localisation of hub genes |
| FinnGen R12 (I9_CAVS_OPERATED) | GWAS summary statistics | MR outcome |
| eQTL Catalogue r7 / GTEx v8 | cis-eQTL | MR instruments |

Raw data are not re-hosted here; download from GEO/FinnGen using the accessions above.

## Repository layout

```
code/
  pipe_common.py            Shared loaders (GEO series matrices, RNA-seq counts), DE, Enrichr, immune scoring
  fetch_gene_lists.py       KEGG immune pathways + UniProt secreted proteins (KW-0964) gene sets
                            (regenerates gene_lists.json used for the intersections)
  analysis_topicB.py        Discovery DE analysis and cross-dataset concordance
  model_topicB.py           LASSO / SVM-RFE / random-forest consensus, diagnostic model, GSE153555 transfer
  full_B.py                 End-to-end pipeline producing the main figures and results_B.json
  upgrade2_mcpcounter.py    MCPcounter deconvolution (official marker set)
  upgrade3_dca.py           Calibration curves, decision-curve analysis, nomogram (CD74+MMP9 model)
  upgrade4_gse12644.py      Evaluation in the fully untouched GSE12644 cohort
  upgrade5_scrna.py         Single-cell localisation of CD74/MMP9 (GSE220774) and stage-wise trends
  upgrade6_stats.py         Reviewer-response statistics: nested cross-validation, bootstrap CIs, calibration slope
  make_fig_upgrade.py       Study workflow figure, DEG heatmap, multi-panel composites
  redo_fig6_B.py            Fig. 6 redraw with DESeq2-based statistics for GSE153555
  run_mr.py / run_mr2.py    Two-sample cis-MR: CD74/MMP9 eQTLs -> FinnGen CAVS
  mr_fetch.py / mr_eqtlc.py / mr_metalcl.py / mr_parquet_probe.py / tabix_remote.py
                            Remote eQTL/GWAS querying helpers (tabix & parquet over HTTP range requests)
  fig7_mr.py                MR forest plot
data/
  mcpcounter_genes.txt      MCPcounter marker genes (ebecht/MCPcounter)
results/
  results_*.json            Numerical outputs cited in the manuscript
  Table_MR_results.csv      MR estimates (Wald ratio / IVW) for CD74 and MMP9
```

## Reproducibility

Python 3.11+ with `numpy`, `pandas`, `scipy`, `scikit-learn`, `matplotlib`. Scripts read raw files from a local `data/` directory (see `pipe_common.py` for expected file names) and write figures/results to an output directory. All randomised steps use fixed seeds.

## Manuscript status

Under preparation for submission. If you use this code, please cite the manuscript and this repository (Zenodo DOI to be added upon release).
