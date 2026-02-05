# MammalianBrainEvolutionAtlas

**Paper Title:** Cross-species single cell atlas of mammalian brain evolution

## Overview
This repository contains the analysis code and scripts for the study "Cross-species single cell atlas of mammalian brain evolution". We integrated single-cell transcriptomic data across multiple mammalian species to construct a comprehensive cell-type atlas. Our analysis investigates cellular homology, species-specific specializations, and the molecular basis of evolutionary phenotypes across the mammalian brain.

## Repository Structure

The code is organized into folders corresponding to the major analysis steps and figures in the manuscript:

### [01_Methods_and_Cortical_Atlas_Construction](./01_Methods_and_Cortical_Atlas_Construction)
Scripts for cross-species data integration, quality control, assembly, and reference mapping.
- **Cross_Species_Integration**: Functions and scripts for integrating human and macaque data.
- **QC_and_Assembly**: Quality control pipelines and genome assembly annotations.
- **Reference_Mapping**: Tools for mapping cells between species (e.g., Macaque to Human).

### [02_L4_Identity_and_Specialization](./02_L4_Identity_and_Specialization)
Analysis focusing on the identity and evolutionary specialization of Layer 4 (L4) cortical neurons.
- Species composition analysis.
- L4 marker gene expression (e.g., *RORB*, *FOXP2*).
- Differential expression analysis between primates and non-primates.

### [03_Thalamic_Atlas_and_Spatial_System](./03_Thalamic_Atlas_and_Spatial_System)
Construction of the cross-species thalamic atlas and investigation of spatial organization.
- Clustering and annotation of thalamic cell types.
- Identification of Midbrain-derived Inhibitory (MDI) neurons.
- Spatial transcriptomics integration (Human and Mouse).
- Correlation of cellular proportions with phenotypic traits (e.g., visual acuity).

### [04_Evolutionary_Phenotypes_and_Genetics](./04_Evolutionary_Phenotypes_and_Genetics)
Downstream analysis linking cellular changes to evolutionary phenotypes and genetics.
- Pathway enrichment analysis.
- Gene expression divergence.
- Correlation with brain weight and cortical folding.

## Citation

If you use this code or data in your research, please cite our paper:
> [Citation Placeholder]

## License

This project is licensed under the MIT License.
