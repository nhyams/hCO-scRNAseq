# Cardiac Organoid scRNAseq Project: Comparing 3 and 4-cell type isogenic cardiac organoids
Workflow of taking 10X Genomics output files: {.features, .barcode, and .matrix} to analyze the impact of adding pericytes into cardiac organoids. This workflow branches out to cover multiple analyses and means to conduct certain analyses (such as preprocessing with seurat vs scanpy, etc).

# Experimental Design
1. iPSC-cardiomyocytes (iCMs) were cultured in 2D, then lifted via TripLE for 40 minutes. Cardiomyocyte suspensions were added into anti-adherant 3D agarose molds. Cardiac spheroids consisting of 3,000 iCMs each formed in agorose molds over the course of 14 days of 3D culture.
2. Supporting cells:
     3-cell type: iPSC-cardiac fibroblasts (iECFs) and iPSC-endothelial cells (iECs)
     4-cell type: iECFs, iECs, and iPSC-pericytes (iPCs)
3. Supporting cells were cultured in 2D, then lifted via TripLE for 10-15 minutes. These were combined in a ratiometric manner on top of cardiac spheroids in a cell suspension.
     3-cell type: 70% iCMs, 15% iECs, 15% iECFs
     4-cell type: 70% iCMs, 15% iECs, 10% iECFs, 5% iPCs
4. Cardiac organoids were allowed to self-aggregate for 5 days of culture.
5. Organoids were dissociated into single cells via TripLE for 40 minutes and assessed for viability and stress markers via flow cytometry.
6. Upon quality-control assessment, 3-cell type and 4-cell type organoids (n = 70 for each group) were processed via the 10X Genomics GEM single cell kit in duplicate (n = 2 libraries for each condition for a total of 4 libraries)-

![Pericytes](https://github.com/user-attachments/assets/06b958a3-7b3e-4da5-b683-fcd89ca43366)
   
# Packages Required

Analysis conducted within R (v4.4.1) were done within the RStudio workspace. 
  R Packages frequently used:
    Seurat
    ggplot2
    tidyverse
    harmony
    CellChat
    NicheNet
    
Analysis conducted within Python (v3.10.7) were done within the VSCode workspace using Jupyter notebooks.
  Python packages used
  Scanpy
  harmonypy
  DecoupleR
  Seaborn
  scvi
