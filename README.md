# Data and Scripts for Statistical Analyses

This repository contains the data and scripts needed to reproduce the statistical analyses conducted for the manuscript:  

*"Shifts in microbial communities driven by the succession of seagrasses to benthic macroalgae may exacerbate heatwave effects in eutrophicated-coastal lagoons"*  
— E. Rubio-Portillo et al.

---

## ⚠️ Note on Heavy Computations

Some of the calculations in this project are very time-consuming. To avoid long compilation times, the corresponding chunks in the `.Rmd` file are set with `eval=FALSE`. Immediately after, there are chunks that **load the precomputed data**.

These data files are located in the `Auxiliary_Data` folder. **Before knitting the `.Rmd` file**, make sure to move them to the **same directory level as the `.Rmd` file**. Notice that the total size of the data files in the `Auxiliary_Data` folder is over 5 GB.

---

## Repository Structure 

* `Auxiliary_Data/` Precomputed data files, to avoid long compilation times
* `Additional-File-2.Rmd` Main RMarkdown file
* `funcionsCODAMETACIRCLE.R` Script with functions used by Additional-File-2.Rmd 
* `Input_Data/` Data obtained in the experiments
* `Main_Results/` Output data tables produced by the statistical analysis
* `README.md` This file

---

## Reproducibility Instructions

1. Clone this repository.
2. Move the contents of `Input_Data` and `Auxiliary_Data` to the same folder as `Additional-File-2.Rmd` and `funcionsCODAMETACIRCLE.R`.
3. Open `Additional-File-2.Rmd` in RStudio and knit.
