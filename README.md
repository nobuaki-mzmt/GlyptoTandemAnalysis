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
remotes::install_version("scales", version = "1.3.0")
remotes::install_version("phytools", version = "2.4-4")
remotes::install_version("maps", version = "3.4.2.1")
remotes::install_version("ape", version = "5.8-1")
remotes::install_version("multcomp", version = "1.4-28")
remotes::install_version("TH.data", version = "1.1-3")
remotes::install_version("MASS", version = "7.3-60.2")
remotes::install_version("mvtnorm", version = "1.3-3")
remotes::install_version("lme4", version = "1.1-36")
remotes::install_version("Matrix", version = "1.7-0")
remotes::install_version("coxme", version = "2.2-22")
remotes::install_version("bdsmatrix", version = "1.3-7")
remotes::install_version("car", version = "3.1-3")
remotes::install_version("carData", version = "3.0-5")
remotes::install_version("survminer", version = "0.5.0")
remotes::install_version("ggpubr", version = "0.6.0")
remotes::install_version("survival", version = "3.6-4")
remotes::install_version("CircMLE", version = "0.3.0")
remotes::install_version("NPCirc", version = "3.1.1")
remotes::install_version("circular", version = "0.5-1")
remotes::install_version("ggridges", version = "0.5.6")
remotes::install_version("viridis", version = "0.6.5")
remotes::install_version("viridisLite", version = "0.4.2")
remotes::install_version("ggplot2", version = "3.5.1")
remotes::install_version("forcats", version = "1.0.0")
remotes::install_version("tidyr", version = "1.3.1")
remotes::install_version("dplyr", version = "1.1.4")
remotes::install_version("data.table", version = "1.17.0")
remotes::install_version("stringr", version = "1.5.1")
remotes::install_version("arrow", version = "19.0.1")
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
