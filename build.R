# Create OmicNavigator study from files in the directory data/

library(openxlsx)
library(dplyr)
library(ggplot2)
library(plotly)
library(OmicNavigator)


# Import data from file

phospho   <- openxlsx::read.xlsx("data/1-s2.0-S0092867420308114-mmc1.xlsx", sheet=1)
protein_abundance <- openxlsx::read.xlsx("data/1-s2.0-S0092867420308114-mmc1.xlsx", sheet=2)

# Create a new study -----------------------------------------------------------

study <- createStudy("SARSCoV2.proteomic.and.phosphoproteomic.profiling",
                     "SARS-CoV-2 global proteomic and phosphoproteomic analysis from Bouhaddou et al converted into a multi-omic OmicNavigator study. In this work Vero E6 cells were profiled at six time points following SARS-CoV-2 infection with quantitative mass spectrometry. Protein abundance and phosphorylation levels were estimated from this data. Please click on the analysis details button for more information.",
                     version = "0.1.0")

# Models -----------------------------------------------------------------------

models <- list(
  protein_abundance = "Differential protein abundance analysis in SARS-CoV-2 infected or mock treated Vero E6 cells",
  protein_phosphorylation = "Differential protein phosphorylation analysis in SARS-CoV-2 infected or mock treated Vero E6 cells" 
)
study <- addModels(study, models)

# Features ---------------------------------------------------------------------

# List of unique uniprot entries in phosphorylation group & having finite nr

cleaned_phospho_features <- which(!duplicated(phospho$C.s.uniprot) &
                                    is.finite(as.numeric(phospho[, 'Ctrl_24Hr.log2FC'])))

cleaned_abundance_features <- which(!duplicated(protein_abundance$C.s.uniprot) &
                                      is.finite(as.numeric(protein_abundance[, 'Ctrl_24Hr.log2FC'])))

# FeatureIDs in the first column are required to be a character vector
# with unique entries
features_phospho <- phospho[, c(4, 1:3, 5, 20:ncol(phospho))]
features_phospho <- features_phospho[cleaned_phospho_features,]
features_phospho <- data.frame(apply(features_phospho, 2, as.character))

features_abundance <- protein_abundance[, c(3, 1:2, 18:ncol(protein_abundance))]
features_abundance <- features_abundance[cleaned_abundance_features,]
features_abundance <- data.frame(apply(features_abundance, 2, as.character))

features <- list(
  protein_abundance = features_abundance[1:3],
  protein_phosphorylation = features_phospho[1:5]
)
study <- addFeatures(study, features)

# MetaFeatures -----------------------------------------------------------------

metaFeatures <- list(
  protein_abundance = features_abundance[, c(1, 4:ncol(features_abundance))],
  protein_phosphorylation = features_phospho[, c(1, 6:ncol(features_phospho))]
)
study <- addMetaFeatures(study, metaFeatures)

# Mapping ----------------------------------------------------------------------

phospho_id <- data.frame("C.s.uniprot" = phospho[cleaned_phospho_features, "C.s.uniprot"])
phospho_id$protein_phosphorylation <- phospho_id$C.s.uniprot

abundance_id <- data.frame("C.s.uniprot" = protein_abundance[cleaned_abundance_features, "C.s.uniprot"])
abundance_id$protein_abundance <- abundance_id$C.s.uniprot

mapping <- dplyr::full_join(x = abundance_id, y = phospho_id, by="C.s.uniprot")[,2:3]

mapping <- list(
  default = mapping
)
study <- addMapping(study, mapping)

# Results (differential expression) --------------------------------------------

phospho_results <- phospho[cleaned_phospho_features, c(4, 6:19)]
phospho_results <- phospho_results %>% mutate_at(c(2:ncol(phospho_results)), as.numeric)

abundance_results <- protein_abundance[cleaned_abundance_features, c(3, 4:17)]
abundance_results <- abundance_results %>% mutate_at(c(2:ncol(abundance_results)), as.numeric)

# The featureID in the first column must be a character vector
results <- list(
  protein_abundance = list(
    ctrl24h_vs_ctrl00h = abundance_results[, c('C.s.uniprot', 'Ctrl_24Hr.log2FC', 'Ctrl_24Hr.adj.pvalue')],
    inf00h_vs_ctrl00h  = abundance_results[, c('C.s.uniprot', 'Inf_00Hr.log2FC', 'Inf_00Hr.adj.pvalue')],
    inf02h_vs_ctrl00h  = abundance_results[, c('C.s.uniprot', 'Inf_02Hr.log2FC', 'Inf_02Hr.adj.pvalue')],
    inf04h_vs_ctrl00h  = abundance_results[, c('C.s.uniprot', 'Inf_04Hr.log2FC', 'Inf_04Hr.adj.pvalue')],
    inf08h_vs_ctrl00h  = abundance_results[, c('C.s.uniprot', 'Inf_08Hr.log2FC', 'Inf_08Hr.adj.pvalue')],
    inf12h_vs_ctrl00h  = abundance_results[, c('C.s.uniprot', 'Inf_12Hr.log2FC', 'Inf_12Hr.adj.pvalue')],
    inf24h_vs_ctrl00h  = abundance_results[, c('C.s.uniprot', 'Inf_24Hr.log2FC', 'Inf_24Hr.adj.pvalue')]
  ),
  protein_phosphorylation = list(
    ctrl24h_vs_ctrl00h = phospho_results[, c('C.s.uniprot', 'Ctrl_24Hr.log2FC', 'Ctrl_24Hr.adj.pvalue')],
    inf00h_vs_ctrl00h  = phospho_results[, c('C.s.uniprot', 'Inf_00Hr.log2FC', 'Inf_00Hr.adj.pvalue')],
    inf02h_vs_ctrl00h  = phospho_results[, c('C.s.uniprot', 'Inf_02Hr.log2FC', 'Inf_02Hr.adj.pvalue')],
    inf04h_vs_ctrl00h  = phospho_results[, c('C.s.uniprot', 'Inf_04Hr.log2FC', 'Inf_04Hr.adj.pvalue')],
    inf08h_vs_ctrl00h  = phospho_results[, c('C.s.uniprot', 'Inf_08Hr.log2FC', 'Inf_08Hr.adj.pvalue')],
    inf12h_vs_ctrl00h  = phospho_results[, c('C.s.uniprot', 'Inf_12Hr.log2FC', 'Inf_12Hr.adj.pvalue')],
    inf24h_vs_ctrl00h  = phospho_results[, c('C.s.uniprot', 'Inf_24Hr.log2FC', 'Inf_24Hr.adj.pvalue')]
  )
)

results$protein_abundance <- lapply(results$protein_abundance, setNames, c('C.s.uniprot', 'log2FC', 'adj.pvalue'))
results$protein_phosphorylation <- lapply(results$protein_phosphorylation, setNames, c('C.s.uniprot', 'log2FC', 'adj.pvalue'))

# sort df by adj pvalue
results$protein_abundance <- lapply(results$protein_abundance, function(x){x[order(x$adj.pvalue),]})
results$protein_phosphorylation <- lapply(results$protein_phosphorylation, function(x){x[order(x$adj.pvalue),]})
  
study <- addResults(study, results)

# Tests ------------------------------------------------------------------------

# here we can add a description of each test, along with columns 
tests <- list(
  protein_abundance = list(
    ctrl24h_vs_ctrl00h = list(
      description = 'Test of differential expression (abundance proteomics) between control 24h and control 0h',
      C.s.uniprot = 'Chlorocebus sabaeus Uniprot ID (and viral protein identifiers)',
      log2FC = 'Log2 transform of the fold change measured as ratio of intensity in Ctrl_24Hr vs Ctrl_00Hr',
      adj.pvalue = 'Adjusted p.value testing the null hypothesis that the log2FC for Ctrl_24Hr vs Ctrl_00Hr is  0.0'
    ),
    inf00h_vs_ctrl00h = list(
      description = 'Test of differential expression (abundance proteomics) between infection 00h and control 0h',
      C.s.uniprot = 'Chlorocebus sabaeus Uniprot ID (and viral protein identifiers)',
      log2FC = 'Log2 transform of the fold change measured as ratio of intensity in Inf_00Hr vs Ctrl_00Hr',
      adj.pvalue = 'Adjusted p.value testing the null hypothesis that the log2FC for Inf_00Hr vs Ctrl_00Hr is  0.0'
    ),
    inf02h_vs_ctrl00h = list(
      description = 'Test of differential expression (abundance proteomics) between infection 02h and control 0h',
      C.s.uniprot = 'Chlorocebus sabaeus Uniprot ID (and viral protein identifiers)',
      log2FC = 'Log2 transform of the fold change measured as ratio of intensity in Inf_02Hr vs Ctrl_00Hr',
      adj.pvalue = 'Adjusted p.value testing the null hypothesis that the log2FC for Inf_02Hr vs Ctrl_00Hr is  0.0'
    ),
    inf04h_vs_ctrl00h = list(
      description = 'Test of differential expression (abundance proteomics) between infection 04h and control 0h',
      C.s.uniprot = 'Chlorocebus sabaeus Uniprot ID (and viral protein identifiers)',
      log2FC = 'Log2 transform of the fold change measured as ratio of intensity in Inf_04Hr vs Ctrl_00Hr',
      adj.pvalue = 'Adjusted p.value testing the null hypothesis that the log2FC for Inf_04Hr vs Ctrl_00Hr is  0.0'
    ),
    inf08h_vs_ctrl00h = list(
      description = 'Test of differential expression (abundance proteomics) between infection 08h and control 0h',
      C.s.uniprot = 'Chlorocebus sabaeus Uniprot ID (and viral protein identifiers)',
      log2FC = 'Log2 transform of the fold change measured as ratio of intensity in Inf_08Hr vs Ctrl_00Hr',
      adj.pvalue = 'Adjusted p.value testing the null hypothesis that the log2FC for Inf_08Hr vs Ctrl_00Hr is  0.0'
    ),
    inf12h_vs_ctrl00h = list(
      description = 'Test of differential expression (abundance proteomics) between infection 12h and control 0h',
      C.s.uniprot = 'Chlorocebus sabaeus Uniprot ID (and viral protein identifiers)',
      log2FC = 'Log2 transform of the fold change measured as ratio of intensity in Inf_12Hr vs Ctrl_00Hr',
      adj.pvalue = 'Adjusted p.value testing the null hypothesis that the log2FC for Inf_12Hr vs Ctrl_00Hr is  0.0'
    ),
    inf24h_vs_ctrl00h = list(
      description = 'Test of differential expression (abundance proteomics) between infection 24h and control 0h',
      C.s.uniprot = 'Chlorocebus sabaeus Uniprot ID (and viral protein identifiers)',
      log2FC = 'Log2 transform of the fold change measured as ratio of intensity in Inf_24Hr vs Ctrl_00Hr',
      adj.pvalue = 'Adjusted p.value testing the null hypothesis that the log2FC for Inf_24Hr vs Ctrl_00Hr is  0.0'
    )
  ),
  protein_phosphorylation = list(
    ctrl24h_vs_ctrl00h = list(
      description = 'Test of differential expression (phosphoproteomics) between control 24h and control 0h',
      C.s.uniprot = 'Chlorocebus sabaeus Uniprot ID (and viral protein identifiers)',
      log2FC = 'Log2 transform of the fold change measured as ratio of intensity in Ctrl_24Hr vs Ctrl_00Hr',
      adj.pvalue = 'Adjusted p.value testing the null hypothesis that the log2FC for Ctrl_24Hr vs Ctrl_00Hr is  0.0'
    ),
    inf00h_vs_ctrl00h = list(
      description = 'Test of differential expression (phosphoproteomics) between infection 00h and control 0h',
      C.s.uniprot = 'Chlorocebus sabaeus Uniprot ID (and viral protein identifiers)',
      log2FC = 'Log2 transform of the fold change measured as ratio of intensity in Inf_00Hr vs Ctrl_00Hr',
      adj.pvalue = 'Adjusted p.value testing the null hypothesis that the log2FC for Inf_00Hr vs Ctrl_00Hr is  0.0'
    ),
    inf02h_vs_ctrl00h = list(
      description = 'Test of differential expression (phosphoproteomics) between infection 02h and control 0h',
      C.s.uniprot = 'Chlorocebus sabaeus Uniprot ID (and viral protein identifiers)',
      log2FC = 'Log2 transform of the fold change measured as ratio of intensity in Inf_02Hr vs Ctrl_00Hr',
      adj.pvalue = 'Adjusted p.value testing the null hypothesis that the log2FC for Inf_02Hr vs Ctrl_00Hr is  0.0'
    ),
    inf04h_vs_ctrl00h = list(
      description = 'Test of differential expression (phosphoproteomics) between infection 04h and control 0h',
      C.s.uniprot = 'Chlorocebus sabaeus Uniprot ID (and viral protein identifiers)',
      log2FC = 'Log2 transform of the fold change measured as ratio of intensity in Inf_04Hr vs Ctrl_00Hr',
      adj.pvalue = 'Adjusted p.value testing the null hypothesis that the log2FC for Inf_04Hr vs Ctrl_00Hr is  0.0'
    ),
    inf08h_vs_ctrl00h = list(
      description = 'Test of differential expression (phosphoproteomics) between infection 08h and control 0h',
      C.s.uniprot = 'Chlorocebus sabaeus Uniprot ID (and viral protein identifiers)',
      log2FC = 'Log2 transform of the fold change measured as ratio of intensity in Inf_08Hr vs Ctrl_00Hr',
      adj.pvalue = 'Adjusted p.value testing the null hypothesis that the log2FC for Inf_08Hr vs Ctrl_00Hr is  0.0'
    ),
    inf12h_vs_ctrl00h = list(
      description = 'Test of differential expression (phosphoproteomics) between infection 12h and control 0h',
      C.s.uniprot = 'Chlorocebus sabaeus Uniprot ID (and viral protein identifiers)',
      log2FC = 'Log2 transform of the fold change measured as ratio of intensity in Inf_12Hr vs Ctrl_00Hr',
      adj.pvalue = 'Adjusted p.value testing the null hypothesis that the log2FC for Inf_12Hr vs Ctrl_00Hr is  0.0'
    ),
    inf24h_vs_ctrl00h = list(
      description = 'Test of differential expression (phosphoproteomics) between infection 24h and control 0h',
      C.s.uniprot = 'Chlorocebus sabaeus Uniprot ID (and viral protein identifiers)',
      log2FC = 'Log2 transform of the fold change measured as ratio of intensity in Inf_24Hr vs Ctrl_00Hr',
      adj.pvalue = 'Adjusted p.value testing the null hypothesis that the log2FC for Inf_24Hr vs Ctrl_00Hr is  0.0'
    )
  )
)
study <- addTests(study, tests)

# Linkouts to external resources for the results table -------------------------
 
resultsLinkouts <- list(
  protein_abundance = list(
    Gene_Name = "https://www.genenames.org/data/gene-symbol-report/#!/symbol/",
    C.s.uniprot = "https://www.uniprot.org/uniprotkb/",
    H.s.uniprot = c("https://www.uniprot.org/uniprotkb/",
                    "https://www.phosphosite.org/uniprotAccAction?id=", 
                    "https://www.proteomicsdb.org/proteomicsdb/#protein/search/query?protein_name=")
  ),
  protein_phosphorylation = list(
    Gene_Name = "https://www.genenames.org/data/gene-symbol-report/#!/symbol/",
    C.s.uniprot = "https://www.uniprot.org/uniprotkb/",
    H.s.uniprot = c("https://www.uniprot.org/uniprotkb/",
                    "https://www.phosphosite.org/uniprotAccAction?id=", 
                    "https://www.proteomicsdb.org/proteomicsdb/#protein/search/query?protein_name=")
  )
)
study <- addResultsLinkouts(study, resultsLinkouts)

# Custom plots -----------------------------------------------------------------

x <- getPlottingData(study, 
                     modelID = "protein_abundance",
                     featureID = c("A0A0D9QUI8"),
                     testID = c("ctrl24h_vs_ctrl00h", "inf00h_vs_ctrl00h",
                                "inf02h_vs_ctrl00h", "inf04h_vs_ctrl00h",
                                "inf08h_vs_ctrl00h", "inf12h_vs_ctrl00h",
                                "inf24h_vs_ctrl00h"))
head(x)

# Singlefeature + Multitest plot (can also be used as multifeature) 
plotMultiTest_singleFeature <- function(x) {
  
  var_x <- data.frame(lapply(x$results, `[`, 2))
  colnames(var_x)<- names(x$results)
  row.names(var_x) <- as.matrix(lapply(x$results, `[`, 1)[[1]])
  var_x <- data.table::data.table(features = rownames(var_x), var_x)
  var_x <- data.table::melt(var_x, measure.vars = c(2:ncol(var_x)))
  
  var_y <- data.frame(lapply(x$results, `[`, 3))
  colnames(var_y) <- names(x$results)
  row.names(var_y) <- as.matrix(lapply(x$results, `[`, 1)[[1]])
  var_y <- data.table::data.table(features = rownames(var_y), var_y)
  var_y <- data.table::melt(var_y, measure.vars = c(2:ncol(var_y)))
  
  df <- merge(var_x, var_y, by=c("variable", "features"))
  colnames(df)[3:4] <- c("log2FC", "adj.pvalue")
  
  cols <- c("darkred", "cadetblue1", "cadetblue2", "cadetblue3", "cadetblue", "cadetblue4", "darkblue")
  cols <- cols[1:length(names(x$results))]
  
  plotly::plot_ly(data = df, x = ~variable, y = ~log2FC,
                  type = "bar",
                  hoverinfo = 'text',
                  text = ~paste('adj.pval: ', round(adj.pvalue, 6), '\nlog2FC: ', log2FC),
                  color = ~variable,
                  colors = cols) %>%
    layout(xaxis = list(title = ''),
           title = unique(df$features),
           showlegend = FALSE)
}
plotMultiTest_singleFeature(x)

# Singlefeature + Multitest + Multimodel plot 
x <- getPlottingData(study,
                     modelID = c(rep("protein_abundance", 7), 
                                 rep("protein_phosphorylation", 7)),
                     featureID = c("A0A0D9QUI8"), 
                     testID = rep(c("ctrl24h_vs_ctrl00h", "inf00h_vs_ctrl00h",
                                "inf02h_vs_ctrl00h", "inf04h_vs_ctrl00h",
                                "inf08h_vs_ctrl00h", "inf12h_vs_ctrl00h",
                                "inf24h_vs_ctrl00h"), 2)) %>% suppressMessages()
head(x)

plotMultiModel_singleFeature <- function(x) {
  
  mm <- x$mapping[complete.cases(x$mapping),]
  models <- names(x)[1:length(names(x))-1]
  
  for (i_models in models) {
    testID <- names(x[[i_models]]$results)
    for (i_test in testID) {
      selected_feats <- which(x[[i_models]]$results[[i_test]][,1] %in% mm[, i_models])
      x[[i_models]]$results[[i_test]] <- x[[i_models]]$results[[i_test]][selected_feats,]
      log2FC <- x[[i_models]]$results[[i_test]]$log2FC
      feat   <- x[[i_models]]$results[[i_test]]$C.s.uniprot
      pval   <- x[[i_models]]$results[[i_test]]$adj.pvalue
      if (!exists("gdf")) {
        gdf <- data.frame(feature = feat, modelID = i_models, testID = i_test, log2FC = log2FC, adj.pvalue = pval)
      } else {
        gdf <- rbind(gdf, data.frame(feature = feat, modelID = i_models, testID = i_test, log2FC = log2FC, adj.pvalue = pval))
      }
    }
  }
  gdf$cols <- c("darkred", "cadetblue1", "cadetblue2", "cadetblue3", "cadetblue", "cadetblue4", "darkblue")
  gdf$groups <- paste(gdf[,"modelID"], gdf[,"testID"], sep = " / ")
  gdf$log2FC <- round(gdf$log2FC, 4)

  p <- ggplot(data=gdf,
              aes(text = paste('adj p.val: ', round(adj.pvalue, 4),
                               '</br> feature: ', feature),
                  x = testID, y = log2FC)) +
    geom_bar(stat='identity', fill=gdf$cols) + 
    facet_grid(modelID~.) +
    theme(axis.text.x = element_text(angle = 90, vjust = .5),
          panel.background = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.y = element_line(size = 0.5, linetype = 'solid',
                                          colour = "gray"), )
  
  ggplotly(p)
}
plotMultiModel_singleFeature(x)

# MultiFeature + singleTest + Multimodel plot

mapped_proteins <- study$mapping$default
mapped_proteins <- mapped_proteins[!is.na(mapped_proteins$protein_abundance) & !is.na(mapped_proteins$protein_phosphorylation), "protein_abundance"]
x <- getPlottingData(study,
                     modelID = c(rep("protein_phosphorylation", 7), 
                                 rep("protein_abundance", 7)),
                     featureID = mapped_proteins,
                     testID = rep(c("ctrl24h_vs_ctrl00h", "inf00h_vs_ctrl00h",
                                    "inf02h_vs_ctrl00h", "inf04h_vs_ctrl00h",
                                    "inf08h_vs_ctrl00h", "inf12h_vs_ctrl00h",
                                    "inf24h_vs_ctrl00h"), 2)) %>% suppressMessages()
head(x)

plotMultiModel_multiFeature <- function(x) {
  
  mm <- x$mapping[complete.cases(x$mapping),]
  if (nrow(mm) == 0) stop('Selected features cannot be mapped across models') 
  
  p_results  <- x$protein_abundance$results[[1]][x$protein_abundance$results[[1]]$C.s.uniprot %in% mm$protein_abundance,]
  t_results  <- x$protein_phosphorylation$results[[1]][x$protein_phosphorylation$results[[1]]$C.s.uniprot %in% mm$protein_phosphorylation,]
  t_features <- x$protein_phosphorylation$features[x$protein_phosphorylation$features$C.s.uniprot %in% mm$protein_phosphorylation,]
  
  # ggdf to store data for plot 
  minuslog10_pval <- -log10(t_results$adj.pvalue)
  threshold       <- t_results$adj.pvalue < 0.05 & abs(t_results$log2FC) > 1.8
  gene_labels     <- as.character(t_features$Gene_Name)
  
  ggdf <- as.data.frame(
    cbind(
      feature     = as.character(gene_labels),
      log2FC_m1   = round(t_results$log2FC, 4),
      log2FC_m2   = round(p_results$log2FC, 4),
      minuslog10_pval = ifelse(minuslog10_pval < 3, 3, minuslog10_pval),
      pval_transc = as.character(ifelse(t_results$adj.pvalue >= 0.001,
                                        round(t_results$adj.pvalue, 3),
                                        "<0.001")),
      pval_proteo = as.character(ifelse(p_results$adj.pvalue >= 0.001,
                                        round(p_results$adj.pvalue, 3),
                                        "<0.001")
      )
    )
  )
  ggdf$log2FC_m1 <- as.numeric(ggdf$log2FC_m1)
  ggdf$log2FC_m2 <- as.numeric(ggdf$log2FC_m2)
  
  plotly::plot_ly(data = ggdf, x = ~log2FC_m1, y = ~log2FC_m2,
                  type = "scatter",
                  mode = "markers",
                  hoverinfo = 'text',
                  text = paste(ggdf$feature, '</br></br>',
                               "phosphorylation adj pval: ", ggdf$pval_transc, '</br>',
                               "phosphorylation log2FC: ", ggdf$log2FC_m1, '</br>',
                               "abundance adj pval: ", ggdf$pval_proteo, '</br>',
                               "abundance log2FC: ", ggdf$log2FC_m2) ,
                  marker = list(size = ~minuslog10_pval)) %>%
    layout(
      yaxis = list(title = paste0('log2FC abundance | ', names(x$protein_abundance$results)[[1]])),
           xaxis = list(title = paste0('log2FC phosphorylation | ', names(x$protein_phosphorylation$results)[[1]])))
}
plotMultiModel_multiFeature(x)

# Cannot have duplicate plots in different models
phosphoplot_MultiTest_singleFeature = plotMultiTest_singleFeature
phosphoplot_MultiModel_singleFeature = plotMultiModel_singleFeature
phosphoplot_MultiModel_multiFeature = plotMultiModel_multiFeature

plots <- list(
  protein_abundance = list(
    plotMultiTest_singleFeature = list(
      displayName = "barplot log2FC/pval per test",
      packages = c("data.table", "plotly"),
      plotType = c("multiTest", "plotly")
    ),
    plotMultiModel_singleFeature = list(
      displayName = "multiModel barplot log2FC/pval per test",
      packages = c("plotly"),
      plotType = c("singleFeature", "multiTest", "multiModel", "plotly")
    ),
    plotMultiModel_multiFeature = list(
      displayName = "multiModel scatterplot",
      packages = c("data.table", "plotly"),
      plotType = c("multiFeature", "singleTest", "multiModel", "plotly")
    )
  ),
  protein_phosphorylation = list(
    phosphoplot_MultiTest_singleFeature = list(
      displayName = "barplot log2FC/pval per test",
      packages = c("data.table", "plotly"),
      plotType = c("multiTest", "plotly")
    ),
    phosphoplot_MultiModel_singleFeature = list(
      displayName = "multiModel barplot log2FC/pval per test",
      packages = c("plotly"),
      plotType = c("singleFeature", "multiTest", "multiModel", "plotly")
    ),
    phosphoplot_MultiModel_multiFeature = list(
      displayName = "multiModel scatterplot",
      packages = c("data.table", "plotly"),
      plotType = c("multiFeature", "singleTest", "multiModel", "plotly")
    )
  )
)

study <- addPlots(study, plots = plots)

plotStudy(study, 
          modelID = "protein_abundance", 
          featureID = "A0A0D9R1T3",
          testID = c("ctrl24h_vs_ctrl00h", "inf00h_vs_ctrl00h",
                     "inf02h_vs_ctrl00h", "inf04h_vs_ctrl00h",
                     "inf08h_vs_ctrl00h", "inf12h_vs_ctrl00h",
                     "inf24h_vs_ctrl00h"),
          plotID = "plotMultiTest_singleFeature") %>% suppressMessages()

plotStudy(study, 
          modelID = c(rep("protein_abundance", 7), rep("protein_phosphorylation", 7)), 
          featureID = "A0A0D9QUI8",
          testID = rep(c("ctrl24h_vs_ctrl00h", "inf00h_vs_ctrl00h",
                         "inf02h_vs_ctrl00h", "inf04h_vs_ctrl00h",
                         "inf08h_vs_ctrl00h", "inf12h_vs_ctrl00h",
                         "inf24h_vs_ctrl00h"), 2),
          plotID = 'plotMultiModel_singleFeature') %>% suppressMessages()
          
plotStudy(study, 
          modelID = c(rep("protein_abundance", 7), rep("protein_phosphorylation", 7)),
          featureID = mapped_proteins,
          testID = rep(c("ctrl24h_vs_ctrl00h", "inf00h_vs_ctrl00h",
                         "inf02h_vs_ctrl00h", "inf04h_vs_ctrl00h",
                         "inf08h_vs_ctrl00h", "inf12h_vs_ctrl00h",
                         "inf24h_vs_ctrl00h"), 2),
          plotID = 'plotMultiModel_multiFeature') %>% suppressMessages()

plotStudy(study, 
          modelID = "protein_phosphorylation", 
          featureID = "A0A0D9R1T3",
          testID = c("ctrl24h_vs_ctrl00h", "inf00h_vs_ctrl00h",
                     "inf02h_vs_ctrl00h", "inf04h_vs_ctrl00h",
                     "inf08h_vs_ctrl00h", "inf12h_vs_ctrl00h",
                     "inf24h_vs_ctrl00h"),
          plotID = "phosphoplot_MultiTest_singleFeature") %>% suppressMessages()

plotStudy(study, 
          modelID = c(rep("protein_phosphorylation", 7), rep("protein_abundance", 7)), 
          featureID = "A0A0D9QUI8",
          testID = rep(c("ctrl24h_vs_ctrl00h", "inf00h_vs_ctrl00h",
                         "inf02h_vs_ctrl00h", "inf04h_vs_ctrl00h",
                         "inf08h_vs_ctrl00h", "inf12h_vs_ctrl00h",
                         "inf24h_vs_ctrl00h"), 2),
          plotID = 'phosphoplot_MultiModel_singleFeature') %>% suppressMessages()

plotStudy(study, 
          modelID = c(rep("protein_phosphorylation", 7), rep("protein_abundance", 7)), 
          featureID = mapped_proteins,
          testID = rep(c("ctrl24h_vs_ctrl00h", "inf00h_vs_ctrl00h",
                         "inf02h_vs_ctrl00h", "inf04h_vs_ctrl00h",
                         "inf08h_vs_ctrl00h", "inf12h_vs_ctrl00h",
                         "inf24h_vs_ctrl00h"), 2),
          plotID = 'phosphoplot_MultiModel_multiFeature') %>% suppressMessages()

# Reports ----------------------------------------------------------------------

reports <- list(
  protein_abundance = "data/report.html",
  protein_phosphorylation = "data/report.html"
)
study <- addReports(study, reports)

# Install study package and start app ------------------------------------------

installStudy(study)
if (interactive()) {
  message("Starting app. Should open in new browser tab.")
  startApp()
}
