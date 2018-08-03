########################################################################
### Phosphos data analysis for David Crottes (EGF + TMEM16A experiment)
########################################################################

############################################################
# Step 1 - Data Clean-up
############################################################
## Set working directory
setwd("C:/Users/Tony Lin/Desktop/Wiita_Lab/Projects/Proteomics_project/Crottes EGF and TMEM16A-dependent signaling pathway/analysis/180802 bR1+2+3/")


## Read in phospho raw file
raw.data = read.delim("data/Phospho (STY)Sites.txt", header = TRUE, stringsAsFactors = FALSE)


## Remove reverse proteins, contaminations, and localization probability < 0.75
require(dplyr)
data = dplyr::filter(raw.data, Reverse != "+") %>%
  dplyr::filter(Potential.contaminant != "+") %>%
  dplyr::filter(Localization.prob >= 0.75) %>%
  dplyr::filter(Delta.score >= 8)


## Calculate log2 of intensity
Intensity.names = grep("^Intensity", names(data), value = TRUE)
LOG.names = sub("^Intensity", "LOG2", Intensity.names)
data[, LOG.names] = lapply(data[, Intensity.names], function(x) log2(x))
rm(Intensity.names, LOG.names)


## Filter columns by names and extract protein IDs
filter_cols = function(df) {
  # df = data frame containing phosphoproteomics data
  require(stringr)
  
  # Extract protein name
  df$UniProtID = df$Protein
  
  # Extract gene name
  df$Gene.Name = gsub(";.*", "", df$Gene.names)
  
  LOG2.names = grep("^LOG2\\..*(?<!___[0-9])$", names(df), value = TRUE, perl = TRUE)
  
  keepCols = c("UniProtID", "Gene.Name", "Position", "Amino.acid", LOG2.names)
  return(df[keepCols])
}
dataC = filter_cols(data)


## Filter rows by valid values
## Group filtering to remove -Inf values (sourced from pulsed_SILAC_ozlem.R)
filter_valids = function(df, conditions, min_count, at_least_one = FALSE) {
  # df = data frame containing LOG2 data for filtering and organized by data type
  # conditions = a character vector dictating the grouping
  # min_count = a numeric vector of the same length as "conditions" indicating the minimum 
  #     number of valid values for each condition for retention
  # at_least_one = TRUE means to keep the row if min_count is met for at least one condition
  #     FALSE means min_count must be met across all conditions for retention
  
  log2.names = grep("^LOG2", names(df), value = TRUE)   # Extract LOG2 column names
  cond.names = lapply(conditions, # Group column names by conditions
                      function(x) grep(x, log2.names, value = TRUE, perl = TRUE))
  
  cond.filter = sapply(1:length(cond.names), function(i) {
    df2 = df[cond.names[[i]]]   # Extract columns of interest
    df2 = as.matrix(df2)   # Cast as matrix for the following command
    sums = rowSums(is.finite(df2)) # count the number of valid values for each condition
    sums >= min_count[i]   # Calculates whether min_count requirement is met
  })
  
  if (at_least_one) {
    df$KEEP = apply(cond.filter, 1, any)
  } else {
    df$KEEP = apply(cond.filter, 1, all)
  }
  
  return(df)  # No rows are omitted, filter rules are listed in the KEEP column
}
dataCF = filter_valids(dataC, 
                       c("shControl_CTRL", "shControl_EGF", "shTMEM_CTRL", "shTMEM_EGF"),
                       min_count = rep(2, 4),
                       at_least_one = TRUE)


## Normalization
## NORMALIZATION METHOD 1: Quantile normalization for each cell line (alternative to global median normalization)
quantile_norm = function(df, conditions, use_keep = TRUE) {
  # df = data frame containing LOG2 columns for normalization
  # conditions = a character vector dictating the grouping
  # use_keep = filter rows using KEEP column prior to normalization
  require(dplyr)
  require(limma)
  
  log2.names = grep("^LOG2", names(df), value = TRUE)   # Extract LOG2 columns
  cond.names = lapply(conditions, # Group column names by conditions
                      function(x) grep(x, log2.names, value = TRUE, perl = TRUE))
  
  if (use_keep) df = df[df$KEEP, ]
  LOG2_df = as.matrix(df[log2.names])
  LOG2_df[!is.finite(LOG2_df)] = NA
  
  for(group in cond.names) LOG2_df[, group] = normalizeQuantiles(LOG2_df[, group])
  
  df[log2.names] = LOG2_df[, log2.names]
  
  return(df)
}
dataCFN = quantile_norm(dataCF,
                        c("shControl_CTRL", "shControl_EGF", "shTMEM_CTRL", "shTMEM_EGF"),
                        use_keep = TRUE)


## Imputation (remove KEEP = FALSE after imputation)
## Impute missing data: Hybrid missing data imputation (sourced from Jae_surfaceome_analysis.R)
hybrid_impute = function(df, conditions, use_keep = TRUE) {
  # df = data frame containing filtered 
  # conditions = a character vector dictating the grouping
  # use_keep = filter rows using KEEP column prior to imputation (WILL AFFECT IMPUTATION OUTCOME IF DATA NOT FILTERED ALREADY)
  require(imputeLCMD)
  require(dplyr)
  
  # Apply KEEP filter
  if (use_keep) df = df[df$KEEP, ]
  
  # Group column names by condition
  log2DF = dplyr::select(df, starts_with("LOG2"))   # Extract LOG2 dataframe
  DF = dplyr::select(df, -starts_with("LOG2"))
  cond.names = lapply(conditions, # Group column names by conditions
                      function(x) grep(x, names(log2DF), value = TRUE, perl = TRUE))
  
  # Create new columns indicating whether the values are imputed
  impute.names = sub("^LOG2", "IMP", unlist(cond.names))
  DF[impute.names] = lapply(unlist(cond.names), function(x) !is.finite(log2DF[[x]]))
  
  # Imputation
  set.seed(1)
  imputeDF = lapply(cond.names,   # Impute each group separately
                    function(x) {
                      tempDF = as.matrix(log2DF[x])
                      modelS = model.Selector(tempDF)
                      impute.MAR.MNAR(tempDF, modelS, method.MAR = "MLE", method.MNAR = "MinProb")
                    })
  imputeDF = do.call(cbind, imputeDF)   # Combine a list of data frames into a single data frame
  
  return(cbind(DF, imputeDF))
}
dataCFNI = hybrid_impute(dataCFN, 
                         c("shControl_CTRL", "shControl_EGF", "shTMEM_CTRL", "shTMEM_EGF"),
                         use_keep = TRUE)


## Perform Welch's t-test (copied from MaxQuant_report_Christine.Rmd)
welch_test = function(df, diff, bioRep = c("1$", "2$", "3$")) {
  # Function only works for pairwise comparisons
  # df = data frame containing log2 intensity data (AFTER EXCLUDING OR IMPUTING MISSING VALUES)
  # diff = find the Ratio by computing element1 - element2
  # bioRep = biological replicate identifiers
  log2.names = grep("^LOG2", names(df), value = TRUE)   #Extract all LOG2. names
  group1 = grep(diff[[1]], log2.names, value = TRUE)
  group2 = grep(diff[[2]], log2.names, value = TRUE)
  group1 = sapply(bioRep, function(x) grep(x, group1, value = TRUE))
  group2 = sapply(bioRep, function(x) grep(x, group2, value = TRUE))
  
  # Perform Welch's two-sided unpaired two-sample t-test
  result = data.frame(Pvalue = numeric(0), MEAN.RATIO = numeric(0))
  for (i in 1:nrow(df)) {
    stats = tryCatch({t.test(unlist(df[i, group1]), unlist(df[i, group2]), alternative = "two.sided",
                            paired = FALSE, var.equal = FALSE)},
                     error = function(e) list(p.value = NA, estimate = NA))
    result[i, ] = c(stats$p.value, -diff(stats$estimate))
  }
  result$LOG.Pvalue = -log10(result$Pvalue)
  
  return(cbind(df, result))
}
CTL = welch_test(dataCFNI, diff = c("shControl_EGF", "shControl_CTRL"))
TMEM = welch_test(dataCFNI, diff = c("shTMEM_EGF", "shTMEM_CTRL"))  
TvC = welch_test(dataCFNI, diff = c("shTMEM_CTRL", "shControl_CTRL"))
TegfvC = welch_test(dataCFNI, diff = c("shTMEM_EGF", "shControl_CTRL"))


## Print tables
write.csv(CTL, "CTL.csv", row.names = F, quote = F)
write.csv(TMEM, "TMEM.csv", row.names = F, quote = F)
write.csv(TvC, "TvC.csv", row.names = F, quote = F)
write.csv(TegfvC, "TegfvC.csv", row.names = F, quote = F)


## Print table for KSEA App
KSEA_table = function(df) {
  # df = output data frame from "welch_test" function after comparing two groups
  data.frame(Protein = df$UniProtID,
             Gene = df$Gene.Name,
             Peptide = rep("NULL", nrow(df)),
             Residue.Both = paste0(df$Amino.acid, df$Position),
             p = rep("NULL", nrow(df)),
             FC = 2^df$MEAN.RATIO)
}
write.csv(KSEA_table(CTL), "KSEA_CTL.csv", row.names = F, quote = F)
write.csv(KSEA_table(TMEM), "KSEA_TMEM.csv", row.names = F, quote = F)    # Need to remove NA in FC
write.csv(KSEA_table(TvC), "KSEA_TvC.csv", row.names = F, quote = F)
write.csv(KSEA_table(TegfvC), "KSEA_TegfvC.csv", row.names = F, quote = F)


## Analyze KSEA tables by Kinase set enrichment analysis
# NetworKIN score cutoff = 2
#shiny::runGitHub('KSEA', 'casecpb')   # Input each table separately (manual)
KSEA_CTL = read.csv("downloadKSEA-CTL.csv", stringsAsFactors = F)
KSEA_TMEM = read.csv("downloadKSEA-TMEM.csv", stringsAsFactors = F)
KSEA_TvC = read.csv("downloadKSEA-TvC.csv", stringsAsFactors = F)
KSEA_TegfvC = read.csv("downloadKSEA-TegfvC.csv", stringsAsFactors = F)



############################################################
# Step 2 - Data Visualization
############################################################
## Histogram of distribution 
hist_impute = function(df, use_keep = TRUE) {
  # df = data frame containing imputed data
  # use_keep = logical indicating whether to use KEEP column to filter rows
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  require(stringr)
  
  if (use_keep) {
    df = filter(df, KEEP)
  }
  
  LOG2.df = dplyr::select(df, starts_with("LOG2"))
  impute.df = dplyr::select(df, starts_with("IMP"))
  
  # Reshape data into columns
  LOG2.df = gather(LOG2.df, "sample", "intensity")
  impute.df = gather(impute.df, "sample", "IMP")
  
  # Combine data
  combine.df = bind_cols(LOG2.df, impute.df["IMP"])
  
  # Create labels
  combine.df = mutate(combine.df, sample = sub("^LOG2\\.", "", sample)) %>%
    mutate(replicate = str_extract(sample, "_bR[1-9]")) %>%
    mutate(replicate = sub("_", "", replicate)) %>%
    mutate(sample = sub("_bR[1-9]*", "", sample))
  
  ggplot(combine.df, aes(x = intensity, fill = IMP)) +
    geom_histogram(alpha = 0.3, binwidth = 0.2, position = "identity") +
    labs(x = expression("log"[2]*"-transformed Intensity"), y = "Frequency") +
    facet_grid(replicate ~ sample) +
    scale_fill_discrete(name = "Imputed",
                        breaks = c("FALSE", "TRUE"),
                        labels = c("-", "+"))
}
#hist_impute(dataCFNI, use_keep = TRUE)


## Pairs Plot
plot_pairs = function(df, pattern, use_keep = FALSE) {
  # df = data frame carrying data of interest
  # pattern = regex pattern to select column of interest
  # use_keep = TRUE means to construct plot on filtered values; FALSE uses all available values
  require(gpairs)
  require(scales)
  if (use_keep) {
    plot.df = df[df$KEEP, grep(pattern, names(df), value = TRUE, perl = TRUE)]
  } else {
    plot.df = df[grep(pattern, names(df), value = TRUE, perl = TRUE)]
  }
  # AUTOMATICALLY EXCLUDE -Inf VALUES
  plot.df = as.matrix(plot.df)
  plot.df[is.infinite(plot.df)] = NA
  colnames(plot.df) = gsub("^LOG2.", "", colnames(plot.df))
  colnames(plot.df) = gsub("_bR", "_", colnames(plot.df))
  
  gpairs(plot.df,
         upper.pars = list(scatter = "lm"),
         scatter.pars = list(pch = 20,
                             col = alpha("black", 0.3)),
         lower.pars = list(scatter = "stats"),
         stat.pars = list(verbose = FALSE, fontsize = 10))
}
#plot_pairs(dataCFN, "^LOG2.shControl_CTRL", TRUE)
#plot_pairs(dataCFN, "^LOG2.shControl_EGF", TRUE)
#plot_pairs(dataCFN, "^LOG2.shTMEM_CTRL", TRUE)
#plot_pairs(dataCFN, "^LOG2.shTMEM_EGF", TRUE)


## Volcano plot (copied from Jae_surfaceome_analysis.R)
plot_volcano = function(df, ratio_cutoff = 1, P_cutoff = 0.05, 
                        use_adjusted = FALSE, labeling = TRUE, title = "") {
  # df = dataframe containing plot data
  # ratio_cutoff = threshold of biological significance
  # P_cutoff = threshold of statistical significance
  # use_adjusted = boolean indicating whether to use Pvalue or FDR-adjusted Pvalue (from pValue_correction function)
  # labeling = boolean indicating whether to include point labels
  # title = character of length 1 indicating plot title
  if (use_adjusted) {
    df$Pvalue = p.adjust(df$Pvalue, method = "BH")
    df$LOG.Pvalue = -log10(df$Pvalue)
  }
  
  df$color = "black"
  df$color[df$MEAN.RATIO > ratio_cutoff & df$Pvalue < P_cutoff] = "darkgreen"
  df$color[df$MEAN.RATIO < -ratio_cutoff & df$Pvalue < P_cutoff] = "red"
  df$trans = 0.1   # point transparency
  df$trans[df$color != "black"] = 0.3
  df$color = factor(df$color)
  
  require(ggplot2)
  require(ggrepel)
  fig = ggplot(df) +
    geom_point(aes(MEAN.RATIO, LOG.Pvalue, colour = color), alpha = df$trans) +
    theme_minimal() +
    scale_colour_manual(values = levels(df$color)) + 
    geom_hline(yintercept = -log10(P_cutoff), linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = c(ratio_cutoff, -ratio_cutoff), linetype = "dashed", color = "gray50") +
    scale_x_continuous(
      name = expression("log"[2]*"( fold change )"),
      breaks = seq(-10, 10, 1)) +
    theme(legend.position="none")
  #annotate("text", x = 1.5, y = -log10(P_cutoff)-0.1, label = paste0("P-value = ", P_cutoff), 
  #         size = 3, colour = "gray50")
  
  if (labeling) {
    fig = fig + geom_text_repel(
      data = subset(df, df$color != "black"),
      aes(MEAN.RATIO, LOG.Pvalue, label = Gene.Name, colour = color),
      size = 3)
  }
  
  if (title != "") {
    fig = fig + ggtitle(title)
  }
  
  if (use_adjusted) {
    print(fig + ylab(expression("- log"[10]*"( adjusted P-value )")))
  } else {
    print(fig + ylab(expression("- log"[10]*"( P-value )")))
  }
}
plot_volcano(CTL, use_adjusted = T, labeling = FALSE)


## HEATMAP VISUALIZATION FUNCTION: pheatmap on Kinase Activities
plot_heatmap = function(mtrx) {
  # mtrx = data matrix to be plotted as heatmap
  require(pheatmap)
  require(dendsort)
  require(RColorBrewer)
  sampleDist = dist(t(mtrx), method = "euclidean")   # Euclidean distances between samples
  sitesDist = dist(mtrx, method = "euclidean")   # Euclidean distances between kinases
  
  sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))  # For ordering dendrogram
  
  #colors = colorRampPalette(rev(brewer.pal(11, "RdBu")))(300)
  colors = colorRampPalette(c("blue", "white", "red"))(300)
  breakPTS = seq(from = -max(abs(mtrx)),
                 to = max(abs(mtrx)),
                 length.out = 300)
  
  pheatmap(t(mtrx),
           cluster_cols = sort_hclust(hclust(sitesDist)),
           cluster_rows = sort_hclust(hclust(sampleDist)),
           col = colors,
           breaks = breakPTS,
           fontsize_row = 10,
           fontsize_col = 10,
           number_color = "black", 
           fontsize_number = 15)
}


## KEGG pathway mapping
## Convert kinase IDs to KEGG ID
require(httr)
# Map uniprotID to KEGG identifiers
GET(url = "http://rest.kegg.jp", path = "conv/uniprot/hsa",
    write_disk("KEGG_uniprot_ID_map.txt", overwrite = TRUE))
IDmap = read.table("KEGG_uniprot_ID_map.txt", sep = "\t", stringsAsFactors = FALSE)
colnames(IDmap) = c("kid", "uniprotID")   # name columns
IDmap$kid = sub("^hsa:", "", IDmap$kid)   # clean characters
IDmap$uniprotID = sub("^up:", "", IDmap$uniprotID)   # clean characters

# Convert kinase names to UniProtID
require(gProfileR)
require(dplyr)
uniprot = gconvert(KSEA_CTL$Kinase.Gene, target = "UNIPROTSWISSPROT",
                   mthreshold = 1, filter_na = F)   # same set of kinase across conditions
uniprot = mutate_if(uniprot, is.factor, as.character) %>%
  dplyr::select(alias, target) %>%
  dplyr::rename(Kinase.Gene = alias, uniprotID = target)


KEGGmapper = function(df, convertUNIPROT, convertKID, pathwayID, title) {
  # df = output data frame from KSEA App
  # convertUNIPROT = data frame containing two columns "Kinase.Gene" and "uniprotID"
  # convertKID = data frame containing two columns "name" and "uniprotID" linking KEGG id to uniprotID
  # pathwayID = a character vector of length 1 indicating KEGG pathway map for visualization
  # title = character of length 1 indicating output plot name
  require(pathview)
  
  df = left_join(df, convertUNIPROT, by = "Kinase.Gene") %>%
    left_join(convertKID, by = "uniprotID")
  score = df$mS
  names(score) = df$kid

  pathview(gene.data = score, pathway.id = pathwayID, 
           species = "hsa", 
           limit = list(gene = max(abs(score)), cpd = 1),
           bins = list(gene = 20, cpd = 20),
           low = list(gene = "blue", cpd = NA),
           mid = list(gene = "gray", cpd = NA),
           high = list(gene = "red", cpd = NA),
           out.suffix = title)
}
KEGGmapper(KSEA_CTL, uniprot, IDmap, "05212", "CTL-pancreatic")
KEGGmapper(KSEA_CTL, uniprot, IDmap, "04012", "CTL-EGFR")
KEGGmapper(KSEA_TMEM, uniprot, IDmap, "05212", "TMEM-pancreatic")
KEGGmapper(KSEA_TMEM, uniprot, IDmap, "04012", "TMEM-EGFR")
KEGGmapper(KSEA_TvC, uniprot, IDmap, "05212", "TvC-pancreatic")
KEGGmapper(KSEA_TvC, uniprot, IDmap, "04012", "TvC-EGFR")
KEGGmapper(KSEA_TegfvC, uniprot, IDmap, "05212", "TegfvC-pancreatic")
KEGGmapper(KSEA_TegfvC, uniprot, IDmap, "04012", "TegfvC-EGFR")


#save.image("crottes.RData")
