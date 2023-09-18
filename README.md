# QuickFixR
Hello world! Hang in there!

Here I present QuickFixR. A small shiny based R packages, which makes R interactive and more user friendly. The package is trying to help researchers working on any **'Omics based data**, which are not super happy about coding at all. Within this packages it is possible to get an overview of data. **...and actually get some great insigths in your data!**

## Introduction to dependencies

The only thing you need is some data in .csv or .txt format. The package needs three types of table to function.

• **Count Table**
    A table which consists of columns as samples and rows as features of any kind, like bacteria (MAGs, AVSs or OTUs), genes, proteins or even metabolites!
    Columns names are sample names and row names are set as feature IDs.

• **Classification Table**   
  A table which consists of columns as layers of feature information (like Kingdom, Phyla, Order, Class, Family, etc for taxonomy of MAGs). Rownames are set as       feature IDs (the same as the rownames in the count table).

• **Sample Information Table**
  I normally say this is the most important thing. Without sample information, fancy data just becomes useless stuff on your harddrive. Well! Columns as variables, like **sample types**, **Location of sampling**, **Instrument**, **batch**, **Tissue types**, **diet**, etc.
  Rownames are sample names (the same as the colnames in the count table).

See examples below.

![alt text](https://github.com/JacobAgerbo/QuickFixR/blob/main/inst/shiny/www/data_example.png)

## Main Features in QuickFixR

This package is mainly combining already established methods from really nice pacakages like _animalcules_, _lmerTest_ and _boral_. I am only trying to provide a more user-friendly approach. These packages should be fully credited for their amazing work, which can be found [**here**](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-021-01013-0), [**here**](https://www.jstatsoft.org/article/view/v082i13), and [**here**](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12514).

**Now to the features!**
• Uploading of your data

• Getting an overview of your data and metadata

• Investigate composition of your data, either raw or by using fancy metrics, like Hill Diversity

• Ordination, including PCA, PCoA and UMAP (fancy!)

• Prevalence testing

• Differential testing

• Correlation with linear mixed effect models (taking into acount random variable)

• Bayesian ordination and regression analysis (quite amazing, must say!)

• Community analysis, using networks

# Installation

This packages requires quite some packages to run, since it is a rather extensive framework. Installing them can be a bugger.
Installation of this package will require R Version >4.1.0.

First you have to download JAGS for bayesian modelling, because BORAL and rjags is depedent on this. Please find a MAC, WINDOWS, and conda solution below.

**MAC**

https://sourceforge.net/projects/mcmc-jags/

**WINDOWS**

https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Windows/

**CONDA**

```{bash}
conda install -c conda-forge jags
```

Furthermore, It can be a good idea to start with this line, fist :)

```{r Installation of dependencies, include = FALSE}   
dependencies <- c("devtools","BiocManager","shiny","shinyjs",
                  "MultiAssayExperiment","ggplot2","plotly",
                  "vegan","dplyr","magrittr","biomformat",
                  "shinythemes","RColorBrewer","decontam",
                  "animalcules","limma", "broom.mixed", "lmerTest", "performance", "gt", "gtExtras","boral","ggboral", "tidyverse", "pbkrtest", "ggiraph", "hilldiv","sva","factoextra", "phyloseq", "ape")


# Install packages not yet installed
#CRAN
installed_packages <- dependencies %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(dependencies[!installed_packages])}
#BiocManager
installed_packages <- dependencies %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
BiocManager::install(dependencies[!installed_packages])}
#Github
installed_packages <- dependencies %in% rownames(installed.packages())
if (installed_packages[21] == FALSE) {
  remotes::install_github("jthomasmock/gtExtras")}
installed_packages <- dependencies %in% rownames(installed.packages())
if (installed_packages[23] == FALSE) {
  remotes::install_github("mbedward/ggboral")}
```
Now if all went smooth, you should be golden!

```{r Installation, include = FALSE}
# Install devtools from CRAN
install.packages("devtools")

# Or the development version from GitHub:
# install.packages("devtools")
devtools::install_github("JacobAgerbo/QuickFixR")
```

## Installation of conda environment and dependencies for QuickFixR

I base this tutorial on conda and therefore miniconda should be installed prior the tutorial, please see link:
*https://docs.conda.io/en/latest/miniconda.html*

First thing we need to do is, creating a conda environment. 

For this you will a config file with all dependencies. This file has already been made and can be downloaded [**here**](https://github.com/EBI-Metagenomics/holofood-course/blob/main/sessions/Metabolomics/QuickFixR.yml). It is called **QuickFixR.yml**.


```
conda env create -f QuickFixR.yml
```

This environment has installed R (>4.1) with several packages, but a few more is needed. 
These packages are not yet to be found on condas channels and therefore we will install them in R

Launch conda environment and subsequently R, by typing:

```
conda activate QuickFixR #activating the environment
R #starting R
```

Now install dependencies

```
dependencies <- c("boral","ggboral", "pbkrtest", "ggiraph", "hilldiv")
installed_packages <- dependencies %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(dependencies[!installed_packages])}
#BiocManager
installed_packages <- dependencies %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
BiocManager::install(dependencies[!installed_packages])}
#Github
installed_packages <- dependencies %in% rownames(installed.packages())
if (installed_packages[2] == FALSE) {
  remotes::install_github("mbedward/ggboral")}
```

Now please install my R package *QuickFixR* 

```
devtools::install_github("JacobAgerbo/QuickFixR")
```

After this you should be golden! And should be able launch the shiny app simply by typing:

```
QuickFixR::QuickFix()
```



## Running Shiny
After installation, which hopefully should be the most difficult part, then you can run a simple commmand to start of the interactive part.
```{r Run QuickFixR, include = FALSE}
QuickFixR::QuickFix()
```

## What is next?

I am planning to incorporate more features in the following nearest future, which includes

• Generalised Linear Models

• Biomarker Prediction

• Batch Correction of data

• Normalisation, like quantile Norm, VSN, etc.

• Clustering, like KNN, PAM, and HCLUST

Please enjoy :)


