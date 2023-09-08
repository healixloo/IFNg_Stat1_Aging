**README: Data Analysis Scripts**

Welcome to the data analysis repository for the manuscript 'IFN-γ-Stat1 Axis Drives Aging-Associated Loss of Intestinal Tissue Homeostasis and Regeneration.'

In this repository, you will find two main code folders, each dedicated to a specific type of data analysis performed in the manuscript.

**1. Bulk RNA Sequencing Analysis (bulkRNAseq):**
   - This folder contains essential scripts for analyzing bulk RNA sequencing data.
   - **RNAseqpipeline R Script:** Used for processing and analyzing bulk RNA sequencing data.
   - **pipeline_RNA-seq Perl Code:** This code is launched by the RNAseqpipeline for specific tasks.
   - **GeneSetBoxPlot R Script:** Utilized for gene set enrichment analysis.
   - **bulkRNAseq_IntestineEpithelial_VivoVitro.R:** Responsible for bulk RNA sequencing analysis of intestine crypts and organoids, contributing to figures 1, 2, S1, and S2 in the manuscript.

**2. Single-Cell RNA Sequencing Analysis (scRNAseq):**
   - In this folder, you'll find R scripts tailored for single-cell RNA sequencing analysis.
   - **scRNAseq_IntestineEpithelial_InVivo.R:** Analyzing single-cell RNA sequencing data from intestine crypts during aging, contributing to figures 3, 7, and S3.
   - **scRNAseq_IntestineEpithelial_Organoids.R:** Analyzing single-cell RNA sequencing data from intestine organoids upon IFN-γ treatment, contributing to figures 3, 5, 7, S4, and S6.
   - **scRNAseq_IntestineImmune_InVivo.R:** Analyzing single-cell RNA sequencing data from intestine resident immune cells during aging, contributing to figures 4 and S5.
   - **Stat1_Identification.R:** Focused on identifying Stat1 regulators in the IFN-γ-driven pathway, contributing to figure 7.

**Data Directories:**
   - The 'data' and 'data_out' folders contain input processed data used in the analysis scripts.
   - BulkRNAseq data: Generated using an in-house RNAseq pipeline detailed in the 'Bulk RNA Sequencing Analysis' folder.
   - scRNAseq data: Processed using The Cell Ranger Software Suite (Version 3.1.0).

**Data Availability:**
   - Raw data for this study has been deposited in the GEO repository under accession numbers GSE129510, GSE129708, GSE169368, GSE129710, GSE174297, and GSE169351.
   - Intermediate data is also available upon request.

Thank you for exploring our data analysis repository. If you have any questions or require further information, please don't hesitate to reach out.
