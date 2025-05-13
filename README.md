# README
## Article Information
This repository provides access to the data and source code used for the manuscript    
**Loss of pair formation predates the evolution of male-less society in termites**  
**Nobuaki Mizumoto, Toshihisa Yashiro, Simon Hellemans**  
**Paper DOI:** [TBA](XXX)

This study investigates the tandem running behavior of three _Glyptotermes_ termite species, including _G. nakajimai_ (both sexual and asexual populations), _G. fuscus_, _G. satsumensis_. The videos were analyzed using the deep-learning posture tracking software, [SLEAP](https://sleap.ai), to quantify tandem running behavior and compare it across species. Also, this study performs phylogenetic comparative analyses of tandem running behavior and mating systems across the whole termite diversity.
This repository includes data and the Python/R scripts.  
Additional data available elsewhere includes, 
- The models for SLEAP are available at TBA.
- GeneBank sequences.

## Table of Contents
This repository includes tracking data, R codes to analyze it, and Python code for video analysis.  
* [README](./README.md)
* [draft](./draft) - draft related files. removed when published.
* [analysis](./analysis)
  * [code](./analysis/code)
    * [data_prep.py](./analysis/code/data_prep.py): Python script to process and clean data from SLEAP tracking outputs.
    * [plot.R](./analysis/code/plot.R): R script for conducting statistical analysis and generating figures.
    * [phylogeny.R](./analysis/code/phylogeny.R): R script for phylogenetic comparative analysis
    * [sleap_model_evaluation.py](./analysis/code/sleap_model_evaluation.py): Python script for sleap model evaluation
  * [data_raw](./analysis/data_raw) - folder containing raw data
  * [data_fmt](./analysis/data_fmt) - folder containing data converted from raw data
  * [output](./analysis/output) - folder containing outputs

## Setup & Dependencies
The scripts of this project is written in R and Python, tested on Windows 11 (64-bit). Following is the environments.

### R Session Info
R version R version 4.4.1 (2024-06-14 ucrt)
Reproducing R Environment
```r
packages <- c(scales="1.3.0", phytools="2.4-4", maps="3.4.2.1", ape="5.8-1", multcomp="1.4-28", TH.data="1.1-3", MASS="7.3-60.2", mvtnorm="1.3-3", lme4="1.1-36", Matrix="1.7-0", coxme="2.2-22", bdsmatrix="1.3-7", car="3.1-3", carData="3.0-5", survminer="0.5.0", ggpubr="0.6.0", survival="3.6-4", CircMLE="0.3.0", NPCirc="3.1.1", circular="0.5-1", ggridges="0.5.6", viridis="0.6.5", viridisLite="0.4.2", ggplot2="3.5.1", forcats="1.0.0", tidyr="1.3.1", dplyr="1.1.4", data.table="1.17.0", stringr="1.5.1", arrow="19.0.1"); for (pkg in names(packages)) remotes::install_version(pkg, version = packages[pkg])
```

### Python Environment
Python 3.11.4
```bash
pip install \
  h5py==3.13.0 \
  numpy==1.25.0 \
  pandas==2.2.3 \
  scipy==1.15.2 \
  feather-format==0.4.1 \
  pillow==11.2.0
```

## Citation
TBA

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact
Nobuaki Mizumoto: nzm0095@auburn.edu
