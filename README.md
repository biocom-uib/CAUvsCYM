# Data and Scripts for Statistical Analyses

This repository contains the data and scripts needed to reproduce the statistical analyses conducted for the manuscript:  

*"Shifts in microbial communities driven by the succession of seagrasses to benthic macroalgae may exacerbate heatwave effects in eutrophicated-coastal lagoons"*  
— E. Rubio-Portillo et al.

---

## ⚠️ Note on Heavy Computations

Some of the calculations in this project are very time-consuming. To avoid long compilation times, the corresponding chunks in the `.Rmd` file are set with `eval=FALSE`. Immediately after, there are chunks that **load the precomputed data**.

These data files are located in the `Auxiliary_Data` folder. **Before knitting the `.Rmd` file**, make sure to move them to the **same directory level as the `.Rmd` file**.

---

## Reproducibility Instructions

1. Clone this repository.
2. Move the contents of `Auxiliary_Data` to the same folder as `Additional-File-2.Rmd`.
3. Open `Additional-File-2.Rmd` in RStudio and knit.
