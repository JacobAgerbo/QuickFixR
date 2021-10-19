# QuickFixR
Hello world! Hang in there! 

Here I present QuickFixR. A small shiny based R packages, which makes R interactive and more user friendly. The package is trying to help researchers working on any Omics-based data, which are not super happy about coding at all. Within this packages it is possible to get an overview of data.

## Introduction to dependencies

The only thing you need is some data in .csv or .txt format. The package needs three types of table to function. 

• **Count Table**
    A table which consists of columns as samples and rows as features of any kind, like bacteria (MAGs, AVSs or OTUs), genes, proteins or even metabolites!
    Columns names are sample names and rownames are set as feature IDs.

• **Classification Table**   
  A table which consists of columns as layers of feature information (like Kingdom, Phyla, Order, Class, Family, etc for taxonomy of MAGs). Rownames are set as       feature IDs (the same as the rownames in the count table). 

• **Sample Information Table**
  I normally say this is the most important thing. Without sample information, fancy data just becomes useless stuff on your harddrive. Well! Columns as variables,   like sample types, Location, Instrument, batch, Tissue types, etc. 
  Rownames are sample names (the same as the colnames in the count table).

See examples below. 

![alt text](https://github.com/JacobAgerbo/QuickFixR/blob/main/inst/shiny/www/data_example.png)

## Main Features in QuickFixR

This package is mainly combining already established methods from really nice pacakages like _animalcules_, _lmerTest_ and _boral_. I am only trying to provide a more user-friendly approach. These packages should be fully credited for their amazing work!

• Uploading of your data

• Getting an overview of your data and metadata

• Investigate composition of your data, either raw or by using fancy metrics, like Hill Diversity

• Ordination, including PCA, PCoA and UMAP (fancy!)

• Prevalence testing 

• Differential testing 

• Correlation with linear mixed effect models (taking into acount random variable)

• Bayesian ordination and regression analysis (quite amazing, must say!)

## Installation

This packages requires quite some packages to run, since it is rather extensive. Installing them can be a bugger.
Installation of this package will require R Version >4.1.0.

```{r Installation, include = FALSE}
# Install devtools from CRAN
install.packages("devtools")

# Or the development version from GitHub:
# install.packages("devtools")
devtools::install_github("JacobAgerbo/QuickFixR")
```

## Running Shiny
After installation, which hopefully should be the most difficult part, then you can run a simple commmand to start of the interactive part.
```{r Run QuickFixR, include = FALSE}
QuickFixR::QuickFix()
```

## What is next?

I am planning to incoorporate more features in the following nearest future, which includes

• Generalised Linear Models

• Biomarker Prediction

• Batch Correction of data

• Normalisation, like quantile Norm, VSN, etc.

• Clustering, like KNN, PAM, and HCLUST

Please enjoy :) 
