## PN microbiome analysis scripts

### This Github page deposit quality trimming, MOTHUR preprocessing, and downstream analysis scripts for manuscript "Nemolizumab treatment normalizes dysregulated microbiome in lesional prurigo nodularis"

quality_trimming.sh: conduct quality trimming on the raw fastq sequences using fastqc to trim bases with quality scores less than 25
mothur.sh: mothur pipeline to process paired trimmed fastq.gz data files, output otu table as well as taxonomy assignments
data_preprocessing.R: process mothur outputs, conduct low abundance taxon removal
descriptive.R: conduct descriptive analysis on genera and phylum in our data, generate figures used for Figure 2
shannon.R: calculate shannon diversity for each sample using rarefication, test for its significance with clinical variables
bray-curtis.R: calculate bray-curtis distance between samples using rarefication, conduct dimension reduction using NMDS, and test for significance between different skin types and batches
Files start with DE: DE_{taxon}_{method}.R: conduct differential abundance/expression analysis
    taxon: taxonomy-level of phylotype unit, genus-/phylum-level are used in this study
    method: method used for differential analysis, ANCOMBC, DESeq2, and limma-voom are used
plot_DE_{genus, phylum}.R: scripts to plot ANCOMBC results and generate Figure 4
mlr_{baseline, nrsw4}.R: machine learning scripts to classify baseline lesional/non-lesion skin types and predict nemolizumab response status using microbiome and demographic data.
