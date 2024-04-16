# I-DIMENSION

I-DIMENSION is a project that integrate Dmg/hgg Methylomic EpigeNetic Spatial transcriptomic prOteomic subtyping System.

More specific, we hypothese that developing a methylo-genomic subtyping system predictive of proteomic compoents of each sample will enable the future stratification of patients into subtypes that are indicative of treatment response. Spatial transcriptomics we can predict how ubiquitous our combination approaches informed by proteomics will be and whether the tumor microenvironment may play a major of current and future treatment strategies.

To validate our hypothesis, we list 3 key aims to achieve the goal:

* Aim1

  * Developing a novel methylo-genomic subtyping system for DMG.
* Aim2

  * Identify the key proteomic signatures (drug targets) of each specimen used in the DMG methylo-genomic subtyping system.
* Aim3

  * Identify the spatial heterogeneity of disease and relate this back to the proteomic signatures.


## Data Collection

We have five modalities for the data resource.

* Tumor genomics data
  * WGS (Whole genome sequencing)
  * String data
* DNA methylation data
  * EPIC array
* Chromatin data
  * ATAC-seq
* RNAdata
  * sc-RNA-seq
* Phospho-and proteomic profiling
  * Proteomics

## Exploratory Data Analysis

EDA can help us explore the data's underlying patterns, relationships, and distributions. 

We can leverage visualization and statistical techiques to gain insights and inform subsequent decision.

For image-based data (i.e Xenium), we can group the image tiles by cell types to understand what the tiles look like, and evaluate the quality of annotation generally, Such as to check whether the nucleus of T cell and B cell tiles are small, fibroblast tiles are thin and long.

For the count-based data (i.e sc-RNA-seq), we can visualize the density of gene expression to observe the distribution (NB or Gaussian distribution).

## Data Cleaning

Raw data is rarely perfect to use. It often contains missing values, outliers, and inconsistencies. To solve these issues, we usually need to leverage the automatic correction techniques such as imputation, outlier detection and quality control. To further improve the quality of dataset, we would mannually correct some annoations based on prior knowledge to ensure the data is consistent and ready for analysis.

## Feature Engeering

Raw features may not be a good or suitable input for the model. How to embed the domain knowledge from expert to select, transform and create features that are relevant to the problem at hand is essential. Good feature engeering can significantly improve the model performance and interpretability.

## Model Building

Starting from some simple models (Auto-Encoder, Multi-Layer Perceptron) to have a general ideas that how capable the model is. And then we can add more useful modules and leverage some tricks to enhance the feature extraction. 

For example, we can introduce Graph Neural 


## Model Evaluation





## Deployment




## Insights and Action
