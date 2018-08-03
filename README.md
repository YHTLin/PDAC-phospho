# PDAC-phospho
Phosphoproteomic analysis of TMEM16A knockdown and EGF treatment on pancreatic cancer line AsPc-1.

## Motivation
The calcium-activated chloride channel (TMEM16A), also known as anoctamin 1 (ANO1), is over-expressed in many cancers. The transmembrane protein is an attractive cancer target owing to its ability to activate many signaling pathways, including the epidermal growth factor receptor (EGFR) signaling, to support cell growth and migration. However, the effect of TMEM16A on cellular phenotype appears to be cell- and tumor-specific (Molecular Cancer 2017). For instance, Wu et al. found that TMEM16A over-expression promotes proliferation in ER-positive, PR-positive, and HER2-negative MCF-7 cells, but inhibits proliferation in ER-negative, PR-negative, and HER2-negative MDA-MB-435S cells (Oncotarget 2017). This suggests that further validation of TMEM16A as a viable target is needed on a per-cell or -tumor basis.

Pancreatic ductal adenocarcinoma (PDAC) is the most common type of pancreatic cancer, and TMEM16A is often over-expressed in PDAC cells. Given that TMEM16A directly interacts with EGFR, which in turn can trigger downstream growth signaling, it was hypothesized that TMEM16A over-expression drives the aggressive cell growth in PDAC. However, a study by Sauter et al. showed that TMEM16A silencing by siRNA did not reduce proliferation in PDAC cell lines (Eur J Physio 2015). In order to dissect the mechanism behind the sustained growth, we aim to probe the phosphoproteome of the PDAC cell line AsPc-1 by mass spectrometry.

## File Description
+ *data folder*: MaxQuant analysis output.
+ *kegg folder*: Figures genereated from mapping kinase activity onto KEGG pathways.
+ *phospho_crottes.R*: R script for data anlysis.
+ *crottes-180802.RData*: Saved workspace from *phospho_crottes.R* and is loaded by *MaxQuant_report_crottes.Rmd*.
+ *MaxQuant_report_crottes.Rmd*: R Markdown file containing figures and descriptions.
+ *MaxQuant_report_crottes.pdf*: PDF generated from *.Rmd* file.
+ *TMEM16A-EGF-phospho-analysis.xlsx*: Output data file from phospho analysis.
