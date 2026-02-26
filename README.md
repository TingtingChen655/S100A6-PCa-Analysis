# Extracellular S100A6 Acts as a Paracrine Mediator Linking Tumor–Bone Crosstalk to Prostate Cancer Skeletal Metastasis

Welcome to the official GitHub repository for our manuscript: **"Extracellular S100A6 acts as a paracrine mediator linking tumor–bone crosstalk to prostate cancer skeletal metastasis"**. 

This repository contains the processed multi-omics datasets (Transcriptomics and Metabolomics) and the custom bioinformatics pipelines used to uncover the molecular mechanisms by which S100A6 regulates prostate cancer (PCa) progression and bone niche remodeling.

---

## 📖 Abstract
Bone metastasis in prostate cancer is driven not only by tumor-intrinsic programs but also by tumor-derived signals that reprogram the bone microenvironment; however, extracellular mediators linking tumor–bone crosstalk to metastatic progression remain incompletely defined. Here, we identify extracellular S100A6 as a tumor-derived paracrine mediator that functionally links tumor signaling with bone niche remodeling during prostate cancer skeletal metastasis. S100A6 expression and secretion were elevated in primary and metastatic prostate cancer tissues, patient serum, and highly metastatic cell lines. Functional studies demonstrated that S100A6 enhances tumor cell migration and metastatic dissemination in vivo, whereas its depletion reduces systemic and skeletal metastatic burden without significantly affecting intrinsic proliferative capacity. Mechanistically, S100A6 was associated with activation of **EGFR–AKT–mTOR signaling** and mesenchymal features in tumor cells. Within the bone microenvironment, extracellular S100A6 promoted osteoclast differentiation while suppressing osteoblast maturation, accompanied by modulation of **Wnt3a/β-catenin signaling**. Consistently, loss of S100A6 attenuated tumor-associated osteolytic bone destruction and preserved trabecular architecture. Collectively, these findings identify extracellular S100A6 as a paracrine mediator linking tumor-intrinsic signaling and dysregulated bone remodeling, highlighting its role in coordinating tumor–bone crosstalk during prostate cancer skeletal metastasis.

---

## 🗂️ Repository Structure

All required data matrices and analysis scripts are organized as follows:

```text
📦 S100A6-PCa-Bone-Metastasis
 ┣ 📂 Data
 ┃ ┣ 📜 RNAseq_processed_matrix.csv        # Normalized transcriptomic read counts (PC3: Control vs. S100A6-KD)
 ┃ ┗ 📜 Metabolomics_processed_matrix.csv  # Processed and normalized metabolite abundance matrix
 ┣ 📂 Scripts
 ┃ ┣ 📜 01_RNAseq_Analysis.R               # Differential expression analysis (DEGs) and pathway enrichment (GO/KEGG)
 ┃ ┣ 📜 02_Metabolomics_Analysis.R         # Differential metabolite screening and metabolic pathway analysis
 ┗ 📜 README.md                            # Project documentation
```

---

## 📊 Data Description

We performed in-house multi-omics sequencing and profiling on human prostate cancer **PC3 cells**. To investigate the role of S100A6, stable knockdown (KD) cell populations were established by transducing PC3 cells with lentiviral vectors encoding short hairpin RNAs (shRNAs) targeting human S100A6, followed by puromycin selection. 

The processed datasets provided in the `Data/` directory include:

1. **`RNAseq_processed_matrix.csv`**: Contains the normalized gene expression matrix comparing Control vs. S100A6-KD PC3 cells. Raw FASTQ files are available upon reasonable request.
2. **`Metabolomics_processed_matrix.csv`**: Contains the relative abundance of annotated extracellular/intracellular metabolites analyzed via LC-MS/MS.

---

## 💻 System Requirements & Dependencies

The bioinformatic analyses and visualizations were performed using **R (version >= 4.1.0)**. 
To run the scripts successfully, please ensure the following R packages are installed:

- `tidyverse` (Data manipulation and visualization)
- `DESeq2` / `edgeR` (Differential gene expression analysis)
- `clusterProfiler` & `org.Hs.eg.db` (GO/KEGG Enrichment analysis)
- `ggplot2` & `pheatmap` (Data visualization)
- `ropls` (for Metabolomics PCA/OPLS-DA analysis)

---

## 🚀 Usage Instructions

1. **Clone this repository** to your local machine:
   ```bash
   git clone https://github.com/[Your-GitHub-Username]/S100A6-PCa-Bone-Metastasis.git
   ```
2. **Set your Working Directory** in R to the root folder of the cloned repository.
3. **Run the scripts** sequentially in the `Scripts/` folder to reproduce the differential analysis, enrichment plots, and multi-omics integration figures presented in the manuscript.

---

## ✉️ Contact

For any questions regarding the code or data, please contact tingtingchen655@gmail.com.