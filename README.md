# QuickFixR
Hello world! Hang in there! 

Here I present QuickFixR. A small shiny based R packages, which makes R interactive and more userfriendly. The package is trying to help researchers working on any Omics-based data, which are not super happy about coding at all. Within this packages it is possible to get an overview of data.

## Introduction to dependencies

The only thing you need is some data in .csv or .txt format. The package needs three types of table to function. 
These includes:
• Count Table
    A table which consists of columns as samples and rows as features of any kind, like bacteria (MAGs, AVSs or OTUs), genes, proteins or even metabolites!
    Columns names are sample names and rownames are set as feature IDs.

• Classification Table    
  A table which consists of columns as layers of feature information (like Kingdom, Phyla, Order, Class, Family, etc for taxonomy of MAGs). Rownames are set as       feature IDs (the same as the rownames in the count table). 

• Sample Information Table
  I normally say this is the most important thing. Without sample information, fancy data just becomes useless stuff on your harddrive. Well! Columns as variables,   like sample types, Location, Instrument, batch, Tissue types, etc. 
  Rownames are sample names (the same as the colnames in the count table).

See examples below. 

![alt text](https://github.com/JacobAgerbo/QuickFixR/blob/main/inst/shiny/www/data_example.png)

## Main Features in QuickFixR
• Uploading of your data
• Getting an overview of your data and metadata
• Investigate composition of your data, either raw or by using fancy metrics, like Hill Diversity
• Ordination, including PCA, PCoA and UMAP (fancy!)
• Testing, including differential testing, correlation with linear mixed effect models (taking into acount random variable), and bayesian ordination and regression   analysis

## What is next?

I am planning to incoorporate more features in the following nearest future, which includes
• Generalised Linear Models
• Biomarker Prediction
• Batch Correction of data
• Normalisation, like quantile Norm, VSN, etc.
• Clustering, like KNN, PAM, and HCLUST

Please enjoy :) 
