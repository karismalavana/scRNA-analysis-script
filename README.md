# Single Cell Analysis of Glioblastoma for the Identification of Therapeutic Targets 🧬🧠

This project focuses on analyzing single-cell RNA sequencing data of glioblastoma to identify genomic markers that could serve as potential therapeutic targets and candidates for immunotherapy. The analysis uncovers cell clusters in glioblastoma and investigates specific gene markers that may influence immune response, offering new insights for targeted treatments.

## Overview

Glioblastoma (GBM) is an aggressive brain tumor known for its heterogeneity, making it difficult to treat. Immunotherapy is a promising approach but has shown limited success with GBM due to its ability to evade immune system responses. By performing single-cell RNA sequencing (scRNA-seq) analysis, this project aims to discover biomarkers that could predict a patient's response to immunotherapy and identify new therapeutic targets.

## Research Question

**What are the cell types comprising a glioblastoma tumor, and can the presence of certain cell types identify potential GBM therapeutic targets and immunotherapy candidates?**

## Key Findings

- **Seven distinct cell clusters** were identified within glioblastoma.
- **Four genomic markers** were selected as potential candidates for immunotherapy treatments:
  - **Cluster 0:** LHFPL3 (Lipoma HMGIC Fusion Partner-Like 3)
  - **Cluster 1:** SPP1 (Secreted Phosphoprotein 1)
  - **Cluster 4:** C3 (Complement C3)
  - **Cluster 5:** SLC38A1 (Solute Carrier Family 38 Member 1)
- Each marker’s role was determined through literature review and functional analysis to assess their potential as therapeutic targets.

## Methodology

The workflow of the project followed these steps:

1. **Data Collection:** The GBM dataset was acquired from the Gene Expression Omnibus (GSE57872).
2. **scRNA-seq Analysis:** Using the R programming language and the Seurat package, the dataset was processed, and cells were clustered based on their gene expression profiles.
3. **Dimensionality Reduction:** PCA (Principal Component Analysis) and UMAP were used to reduce the dimensionality of the dataset, allowing for visualization of gene clusters.
4. **Marker Identification:** Specific gene markers were identified for each cluster based on their expression levels and fold changes.
5. **Cluster Functional Analysis:** The function of each marker was researched to understand its role in GBM and its potential in immunotherapy.

## Tools and Technologies

- **R**: Statistical computing and graphics
- **Seurat**: An R package for single-cell RNA-seq data analysis
- **Gene Expression Omnibus (GEO)**: A database for gene expression data
- **UMAP and PCA**: Dimensionality reduction techniques for cluster visualization

## Results

The scRNA-seq analysis identified seven clusters in the GBM dataset. Four of these clusters, associated with the LHFPL3, SPP1, C3, and SLC38A1 markers, were identified as potential immunotherapy targets. Further studies are required to confirm these findings and explore the potential therapeutic applications.

## Future Directions

- **Functional assays** to validate the role of identified markers in immunotherapy response.
- **Further analysis** using new GBM datasets to confirm the results.
- **Clinical trials** to investigate the therapeutic potential of these markers in immunotherapy for GBM.

## Contributing

We welcome contributions to improve the analysis, explore additional datasets, or propose new therapeutic candidates. If you're interested in contributing, feel free to submit a pull request or open an issue.

## Acknowledgements

This project was made possible by the guidance of Dr. Abrar Choudhury and the resources provided by the National Center for Biotechnology Information. Special thanks to the research community for advancing the field of cancer genomics and immunotherapy.
