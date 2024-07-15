# Intro
Differential methylation analysis is powered by two comprehensive pipeline, `RnBeads` and customised `Cross-packages` pipeline.
The `RnBeads` is used to generate a quick general overview of methylation data. We developed a `Cross-packages` analytical approaches to perform deeper analysis.

# RnBeads workflow
RnBeads includes six components to conduct the comprehensive analysis.
* [Data import]()
* [Quality control]()
* [Preprocessing]()
* [Covariate inference]()
* [Exploratory analysis]()
* [Differential methylation]()


# Cross-packages workflow
We develop a `cross-packages` analytical approach which is tailered for `EPICv2` array. The outlline of workflow is listed below:

* Data loading
* Quality control
* Preprocessing
  * Pre-filtering
  * Normalisation
  * Post-filtering
  * SVD analysis
  * Batch correction
* Probe-wise differential analysis
  * Differential CpG sites identification
  * Sample-wise methylation visualisation
  * Probe-wise volcano plot
  * Heatmap of differential methylated CpG sites
  * Pathway enrichment analysis
* Differential methylation analysis of regions
* Differential variability analysis
  * Differential CpG sites identification
  * Sample-wise methylation visualisation
  * Probe-wise volcano plot
  * Heatmap of differential methylated CpG sites
  * Pathway enrichment analysis

We have done some differential analysis including [DMG vs Non DMG](./DMG_vs_Non_DMG.html), [DMG vs Non DMG brain tumor](./DMG_vs_Non_DMG_BrainTumor.html), [DMG vs Normal](./DMG_vs_Normal.html), [DMG vs Giolblastoma](./DMG_vs_Giolblastoma.html).



