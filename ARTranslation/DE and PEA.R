#' ---
#' title: "Differential expression analysis and pathway enrichment analysis of LNCaP cell line treatments"
#' ---

#'## Loading packages

library(BiocManager) 
library(ggplot2)     
library(tidyr)       
library(dplyr)
library(ggrepel)
library(openxlsx)

#'## Paths

# create shortcut paths
input <- "input_data/"
output_modified_data <- "output_data/0_modified_data/"
output_imputation <- "output_data/1_Imputation/"
output_DE_results <- "output_data/2_DE_results/"
output_Volcano_Plot <- "output_data/3_Volcano_Plots/"
output_GSEA_Reactome <- "output_data/4_GSEA_Reactome_results/"

#'## Import data

#' The Data contains values that are not a number (NaN) and blanks.
#' They are now assigned "NA" values.

raw_data <- read.delim(
  file = paste0(input, "input_data.txt"), 
  comment.char = "#",
  na.strings = c("", "NaN")
)

#'## Modify data

#' Since single cells contain more than one gene symbol (isoforms,...) 
#' separated by semicolons, we use separate_rows from the tidyr package 
#' to solve this issue. Then we filter for NA vales within Gene symbols and
#' check for distinct Genes within the Genes column.

data_filtered <- raw_data |>
  tidyr::separate_rows(Genes, 
                       sep = ";") |>
  dplyr::select(1:16, Genes) |>
  dplyr::filter(!is.na(Genes))

# filter for unique Gene names, keeps the first occurence
data_filtered <- data_filtered |>
  dplyr::distinct(Genes, 
           .keep_all = TRUE)

# save the modified data in an .xlsx file.
wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb = wb, 
                       sheetName = "modified_data")

openxlsx::writeData(wb = wb,
                    sheet = paste0("modified_data"), 
                    x = data_filtered)

openxlsx::saveWorkbook(wb = wb, 
                       file = paste0(output_modified_data, "modified_data.xlsx"),
                       overwrite = TRUE)

#'## Imputation

#' Now we impute missing values using missForest.
#' Before the imputation takes place, we first filter for NA values:
#' max 1 NA in each condition.

# Function to count NAs in each group and the whole row
count_nas <- function(row) {
  row <- as.numeric(row)
  CTRL_NAs <- sum(is.na(row[1:4]))
  R1881_NAs <- sum(is.na(row[5:8]))
  Enza_NAs <- sum(is.na(row[9:12]))
  C6_NAs <- sum(is.na(row[13:16]))
  Total_NAs <- sum(is.na(row))
  
  return(list(CTRL_NAs = CTRL_NAs, 
              R1881_NAs = R1881_NAs,
              Enza_NAs = Enza_NAs, 
              C6_NAs = C6_NAs, 
              Total_NAs = Total_NAs))
}

# Apply the function to each row and filter based on the conditions
data_for_imputation <- data_filtered |>
  dplyr::rowwise() |>
  dplyr::mutate(NA_counts = list(count_nas(c_across(CTRL_1:C6_4)))) |>
  dplyr::filter(
    NA_counts$CTRL_NAs <= 1 &
      NA_counts$R1881_NAs <= 1 &
      NA_counts$Enza_NAs <= 1 &
      NA_counts$C6_NAs <= 1 &
      NA_counts$Total_NAs <= 4
  ) |>
  dplyr::select(-NA_counts)

# select columns containing only numbers
data_numbers <- as.matrix(data_for_imputation[, 1:16])

rownames(data_numbers) <- data_for_imputation$Genes

# Imputation

# Open a file connection for writing messages
message_file <- file(paste0(output_imputation, "Imputation.txt"), 
                     open = "wt")

sink(file = message_file, 
     append = TRUE, 
     type = "output")

sink(file = message_file, 
     append = TRUE, 
     type = "message")

library(missForest)

# starting time
current_timezone <- Sys.timezone()

starting_time <- format(x = Sys.time(), 
                        tz = current_timezone, 
                        format = "%Y-%m-%d %H:%M:%S")

print(starting_time)

# print begin message
cat("starting Imputation \n")

# imputation
imputed_missForest <- missForest::missForest(xmis = data_numbers,
                                             verbose = TRUE)$ximp

# gene information to imputation

imputation_data <- as.data.frame(imputed_missForest)

imputation_data$Genes <- rownames(imputation_data)

# saving results in an excel file => this doesn't seem to work
wb_Imputation <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb = wb_Imputation, 
                       sheetName = "Imputation_Output")

openxlsx::writeData(wb = wb_Imputation, 
                    sheet = "Imputation_Output", 
                    x = imputation_data)

openxlsx::saveWorkbook(wb = wb_Imputation, 
                       file = paste0(output_imputation, "Imputation.xlsx"),
                       overwrite = TRUE)

# saving results .RData
save(imputation_data, 
     file = paste0(output_imputation, "Imputation.Rdata"))

cat("finished Imputation \n")

# finishing time
end_time <- format(x = Sys.time(), 
                   tz = current_timezone, 
                   format = "%Y-%m-%d %H:%M:%S")

print(end_time)

# stop redirecting messages
sink(type = "output")
sink(type = "message")

# Close the file connection
close(message_file)  


#'## Differential Expression Analysis with limma

library(limma)

# create design matrix
groups <- factor(rep(c("CTRL", "R1881", "Enza", "C6"), each = 4), 
                 levels = c("CTRL", "R1881", "Enza", "C6" ))

design <- model.matrix(~0 + groups)

colnames(design) <- levels(groups)

# define weights
aw <- limma::arrayWeights(object = imputed_missForest, 
                          design = design)

# fit linear model
fit <- limma::lmFit(object = imputed_missForest, 
                    design = design, 
                    weights = aw)

# design contrasts
contrast_matrix <- limma::makeContrasts(
  "R1881 vs CTRL" = R1881 - CTRL,
  "Enza vs CTRL" = Enza - CTRL,
  "C6 vs CTRL" = C6 - CTRL,
  levels = design
  )

fit <- limma::contrasts.fit(fit = fit, 
                            contrasts = contrast_matrix)

# eBayes
fit <- limma::eBayes(fit = fit)

# saving results in a list containing the contrasts
results_list <- list()
for (g in dimnames(fit$contrasts)$Contrasts) {
  results <- limma::topTable(fit = fit, 
                             number=Inf, 
                             coef = g, 
                             adjust.method = "BH")
  
  results$Genes <- rownames(results)
  rownames(results) <- NULL
  results_list[[g]] <- results
}

# saving results in an excel file
wb_DE <- openxlsx::createWorkbook()

for (o in names(results_list)) {
  openxlsx::addWorksheet(wb = wb_DE, 
                         sheetName = o)
  
  openxlsx::writeData(wb = wb_DE, 
                      sheet = o, 
                      x = results_list[[o]])
}
openxlsx::saveWorkbook(wb = wb_DE, 
                       file = paste0(output_DE_results, "DE_output.xlsx"),
                       overwrite = TRUE)

#'## Volcano Plots of DE Analysis

#' Visualization of differentially expressed genes

# create column with gene type up, down, ns
for (i in 1:length(results_list)) {
  results_list[[i]] <- results_list[[i]] |>
    dplyr::mutate(gene_type = dplyr::case_when(logFC >= 0 & adj.P.Val <= 0.05 ~ "up",
                                               logFC <= 0 & adj.P.Val <= 0.05 ~ "down",
                                               TRUE ~ "not significant"))
}

# set colors
cols <- c("up" = "#ffad73", 
          "down" = "#26b3ff", 
          "not significant" = "grey") 

# label points of interest - Androgen Receptor (AR)
AR_list <- list()
for (i in names(results_list)) {
  AR <- results_list[[i]] |>
    filter(Genes %in% c("AR"))
  AR_list[[i]] <- AR
}

# Volcano plot
for (b in names(results_list)) {
  volcanoplot <- ggplot2::ggplot(data = results_list[[b]],
                                 mapping = aes(x = logFC, 
                                               y = -log10(adj.P.Val))) + 
    ggplot2::geom_point(aes(colour = gene_type), 
                        shape = 16,
                        size = 1.5) +
    ggplot2::geom_point(data = AR_list[[b]],
                        shape = 21,
                        size = 2, 
                        fill = "firebrick", 
                        colour = "black") +
    ggplot2::geom_hline(yintercept = -log10(0.05),
                        linetype = "dashed") +
    ggplot2::geom_vline(xintercept = log2(1),
                        linetype = "dashed") +
    ggrepel::geom_label_repel(data = AR_list[[b]],
                              aes(label = Genes),
                              force = TRUE,
                              nudge_y = 0,
                              nudge_x = 0.25) +
    ggplot2::scale_x_continuous(breaks = c(seq(-6, 6, 1)),
                                limits = c(-6, 6)) +
    ggplot2::scale_y_continuous(breaks = c(seq(0, 6, 1)),
                                limits = c(0, 6)) +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::labs(title = gsub(pattern = "_",
                               replacement = " ",
                               x = paste0("Volcano Plot: ", b)),
                  x = "log2(fold change)",
                  y = "-log10(adjusted p-value)",
                  colour = "Expression \nchange") +
    ggplot2::theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
                   legend.title = element_text(size = 12, 
                                               hjust = 0.5),
                   legend.text = element_text(size = 10),
                   plot.title = element_text(size = 13,
                                             hjust = 0.5),
                   axis.title.y = element_text(size = 10),
                   axis.title.x = element_text(size = 10),
                   text = element_text(family = "sans"),
                   panel.background = element_blank())
  
  ggplot2::ggsave(filename = paste0(output_Volcano_Plot, "VP ", b, ".tiff"), 
                  plot = volcanoplot,
                  width = 5.5,
                  height = 6,
                  dpi = 800)
}

#'##  Gene Set Enrichment Analysis Reactome

# Open a file connection for writing messages

message_file <- file(paste0(output_GSEA_Reactome, "GSEA Reactome.txt"), 
                     open = "wt")

sink(message_file, 
     append = TRUE, 
     type = "output")

sink(message_file, 
     append = TRUE, 
     type = "message")

# starting time
current_timezone <- Sys.timezone()

starting_time <- format(x = Sys.time(), 
                        tz = current_timezone, 
                        format = "%Y-%m-%d %H:%M:%S")

print(starting_time)

library(ReactomePA)
library(clusterProfiler)

GSEA_Reactome_list <- list()

for (j in names(results_list)) {
  # create workbooks for saving data in excel
  wb <- openxlsx::createWorkbook()
  
  # print begin message
  cat("begin ", j)
  
  # beginning  time
  time <- format(x = Sys.time(), 
                 tz = current_timezone, 
                 format = "%Y-%m-%d %H:%M:%S")
  
  print(time)
  
  res <- results_list[[j]]
  
  # transform Symbols into Entrez-IDs, which is required for gsePathway
  ids <- clusterProfiler::bitr(res$Genes, 
                               fromType ="SYMBOL", 
                               toType = "ENTREZID", 
                               OrgDb = "org.Hs.eg.db", 
                               drop = FALSE)
  
  # attach Entrez-IDs to res
  res$Entrez <- ids[,"ENTREZID"]
  
  # remove NA
  res <- res[complete.cases(res), ]
  
  # ranking: (logFC * adjusted p value) with a sign according to pos. or neg. logFC
  rankings <- sign(res$logFC)*(-log10(res$adj.P.Val))
  
  # genes as names
  names(rankings) <- res$Entrez
  
  # sort in decreasing order
  rankings <- sort(rankings, decreasing = TRUE)

  GSEA_Reactome <- gsePathway(geneList = rankings,
                              organism = "human",
                              exponent = 1,
                              minGSSize = 10,
                              maxGSSize = Inf,
                              eps = 1e-300,
                              pvalueCutoff = 1,
                              pAdjustMethod = "BH",
                              verbose = TRUE,
                              nPermSimple = 10000,
                              by = "fgsea")
  
  print(warnings())
  
  save(GSEA_Reactome, 
       file = paste0(output_GSEA_Reactome, "GSEA Reactome ", j, ".Rdata"))
  
  # save output in a list for visualization
  GSEA_Reactome_list[[j]] <- GSEA_Reactome
  
  # add worksheets
  openxlsx::addWorksheet(wb = wb, 
                         sheetName = j)
  
  openxlsx::writeData(wb = wb,
                      sheet = j, 
                      x = GSEA_Reactome)
  
  # print how far the loop got => name of DE
  cat("finished ", j)
  
  # finishing time
  current_time <- Sys.time()
  time <- format(current_time, tz = current_timezone, format = "%Y-%m-%d %H:%M:%S")
  print(time)
  
  # saving the workbook
  openxlsx::saveWorkbook(wb = wb, 
                         file = paste0(output_GSEA_Reactome, "GSEA Reactome ", j,".xlsx"),
                         overwrite = TRUE)
}

# finishing time for all DE
end_time <- format(x = Sys.time(), 
                   tz = current_timezone, 
                   format = "%Y-%m-%d %H:%M:%S")

print(end_time)
  
# print the starting time and finishing time as well as the difference to see how long the process took
time_difference <- as.POSIXct(end_time, tz = current_timezone) - as.POSIXct(starting_time, tz = current_timezone)

print(paste("startet", starting_time, "and finished", end_time))

print(time_difference)

# stop redirecting messages
sink(type = "output")
sink(type = "message")

# Close the file connection
close(message_file)  

#' Visualization of pathways

# filter significant pathways
R <- GSEA_Reactome_list[["R1881 vs CTRL"]]@result |>
  dplyr::filter(p.adjust < 0.05)

E <- GSEA_Reactome_list[["Enza vs CTRL"]]@result |>
  dplyr::filter(p.adjust < 0.05)


# get top level pathways from Reactome database
get_top_level_list <- list(R, E)

names(get_top_level_list) <- c("R", "E")

library(ReactomeContentService4R)

for (j in names(get_top_level_list)) {
  for (i in 1:length(get_top_level_list[[j]]$ID)) {
    tryCatch(
      expr = {
        get_top_level_list[[j]][i, "top_level"] <- dplyr::select(ReactomeContentService4R::getPathways(id =  get_top_level_list[[j]]$ID[i],
                                                                                                       species = "Homo sapiens",
                                                                                                       top.level = TRUE),
                                                                 displayName)
      },
      error = function(e) {
        print(e)
        get_top_level_list[[j]][i, "top_level"] <- "NA"
      },
      warning = function (w) {
        print(w)
        get_top_level_list[[j]][i, "top_level"] <- "NA"
      },
      finally = {
        message(paste0("[", i, "] ID processed: ",  get_top_level_list[[j]]$Description[i]), " ",  get_top_level_list[[j]]$ID[i])
      }
    )
  }
}


merged_R_and_E <- dplyr::full_join(x = get_top_level_list[[1]],
                                   y = get_top_level_list[[2]], 
                                   by = "ID",
                                   suffix = c("_RvsC", "_EvsC")) |>
  dplyr::select(ID, 
                Description_RvsC,
                NES_RvsC, 
                top_level_RvsC, 
                Description_EvsC,
                NES_EvsC,
                top_level_EvsC) |>
  dplyr::rename(RvsC = NES_RvsC,
                EvsC = NES_EvsC) |>
  dplyr::mutate(top_level = dplyr::coalesce(top_level_RvsC, top_level_EvsC)) |>
  dplyr::select(-top_level_RvsC, -top_level_EvsC) |>
  dplyr::mutate(Description = dplyr::coalesce(Description_RvsC, Description_EvsC)) |>
  dplyr::select(-Description_RvsC, -Description_EvsC)

# pathway heatmap
library(ComplexHeatmap)
library(circlize)

heatmap_data <- as.matrix(merged_R_and_E[, c("RvsC", "EvsC")])
colnames(heatmap_data) <- c("1:R1881 vs Control", "2:Enza vs Control")
rownames(heatmap_data) <- merged_R_and_E$Description

# Define color function
col_fun <- circlize::colorRamp2(breaks = c(-3, 0, 3), 
                                colors = c("dodgerblue3", "white", "firebrick3"))

# Create the heatmap
png(paste0(output_GSEA_Reactome, "pathway heatmap.tiff"), 
    width = 15,
    height = 25,
    units = "cm", 
    res = 1000)

ht <- ComplexHeatmap::Heatmap(
  heatmap_data, 
  name = "Normalized Enrichment Score",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,   
  show_column_names = FALSE,
  
  show_heatmap_legend = TRUE,
  
  row_split = unlist(merged_R_and_E$top_level),    # Split rows by top_level pathway
  column_split = colnames(heatmap_data),
  
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_title_rot = 0,
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  
  column_names_side = "top",
  
  row_gap = unit(4, "mm"),
  column_gap = unit(3, "mm"),
  
  na_col = "grey",
  
  heatmap_legend_param = list(direction = "horizontal",
                              title_position = "topleft",
                              labels_gp = gpar(fontsize = 10),
                              title_gp = gpar(fontsize = 10, 
                                              fontface = "bold"),
                              legend_width = unit(2.8, "cm"))
)
na_legend <- ComplexHeatmap::Legend(at = 1,
                                    labels = "NA values",
                                    legend_gp = gpar(fill = "grey"))

ht <- ComplexHeatmap::draw(ht, 
                         heatmap_legend_side = "bottom")

ht <- ComplexHeatmap::draw(na_legend,
                         x = unit(11.9, "cm"), 
                         y = unit(0.74, "cm"))

dev.off()

# table with translation pathways
library(gt)
library(gtExtras)

Events_Hierarchy <- ReactomeContentService4R::getEventsHierarchy("Homo sapiens")

translation <- c(Events_Hierarchy$children[[19]][1,]$name,
                 
                 Events_Hierarchy$children[[19]]$children[[1]]$name,
                 
                 Events_Hierarchy$children[[19]]$children[[1]]$children[[1]]$name,
                 Events_Hierarchy$children[[19]]$children[[1]]$children[[1]]$children[[1]]$name,
                 Events_Hierarchy$children[[19]]$children[[1]]$children[[1]]$children[[2]]$name,
                 
                 Events_Hierarchy$children[[19]]$children[[1]]$children[[2]]$name,
                 Events_Hierarchy$children[[19]]$children[[1]]$children[[2]]$children[[1]]$name,
                 Events_Hierarchy$children[[19]]$children[[1]]$children[[2]]$children[[1]]$children[[1]]$name,
                 Events_Hierarchy$children[[19]]$children[[1]]$children[[2]]$children[[1]]$children[[2]]$name,
                 Events_Hierarchy$children[[19]]$children[[1]]$children[[2]]$children[[1]]$children[[3]]$name,
                 Events_Hierarchy$children[[19]]$children[[1]]$children[[2]]$children[[1]]$children[[3]]$children[[5]]$name,
                 Events_Hierarchy$children[[19]]$children[[1]]$children[[2]]$children[[1]]$children[[4]]$name,
                 Events_Hierarchy$children[[19]]$children[[1]]$children[[2]]$children[[1]]$children[[5]]$name,
                 Events_Hierarchy$children[[19]]$children[[1]]$children[[2]]$children[[1]]$children[[6]]$name,
                 
                 Events_Hierarchy$children[[19]]$children[[1]]$children[[2]]$children[[2]]$name,
                 
                 
                 Events_Hierarchy$children[[19]]$children[[1]]$children[[3]]$name,
                 
                 Events_Hierarchy$children[[19]]$children[[1]]$children[[4]]$name,
                 Events_Hierarchy$children[[19]]$children[[1]]$children[[4]]$children[[3]]$name,
                 
                 Events_Hierarchy$children[[19]]$children[[1]]$children[[5]]$name,
                 
                 Events_Hierarchy$children[[19]]$children[[1]]$children[[6]]$name,
                 Events_Hierarchy$children[[19]]$children[[1]]$children[[6]]$children[[1]]$name,
                 Events_Hierarchy$children[[19]]$children[[1]]$children[[6]]$children[[2]]$name,
                 Events_Hierarchy$children[[19]]$children[[1]]$children[[6]]$children[[3]]$name)

library(svglite)
library(webshot2)

translation_table <- merged_R_and_E[,c("ID", "RvsC", "EvsC", "Description")] |>
  dplyr::filter(merged_R_and_E$Description %in% translation) |>
  dplyr::select(ID, Description, RvsC, EvsC) |>
  dplyr::mutate(RvsC = round(RvsC, 2))

translation_table <- translation_table |>
  gt::gt() |>
  gt::tab_header(title = "NES values of Translation pathways") |>
  gtExtras::gt_plt_bar(column = RvsC,
                       color = "lightgoldenrod2",
                       text_color = "black",
                       scale_type = "number",
                       accuracy = 0.01)

# Save the gt table
gt::gtsave(translation_table, 
           filename = "translation table.png",
           path = output_GSEA_Reactome,
           vwidth = 10000,
           vheight = 10000,
           zoom = 10)

#'## session and packages
#' Here is important information about packages and session

# Open a file connection for writing messages
message_file <- file("session and packages.txt", 
                     open = "wt")

sink(file = message_file, 
     append = TRUE, 
     type = "output")

sink(file = message_file, 
     append = TRUE, 
     type = "message")

# starting time
current_timezone <- Sys.timezone()

starting_time <- format(x = Sys.time(), 
                        tz = current_timezone, 
                        format = "%Y-%m-%d %H:%M:%S")
print(starting_time)

sessionInfo()

R <- RStudio.Version()
print("RStudio:")
print(R)

library(purrr)
c("BiocManager", 
  "ggplot2", 
  "tidyr",
  "ggrepel",
  "dplyr",
  "openxlsx",
  "missForest", 
  "limma",
  "org.Hs.eg.db", 
  "clusterProfiler", 
  "ReactomePA",
  "DOSE",
  "ReactomeContentService4R",
  "ComplexHeatmap",
  "circlize",
  "gt",
  "gtExtras",
  "svglite",
  "webshot2",
  "purrr"
  ) |>
  purrr::map(citation) |>
  print(style = "text")

java_version <- system("java -version", 
                       intern = TRUE)

print(java_version, 
      sep = "\n")

# stop redirecting messages
sink(type = "output")
sink(type = "message")

# Close the file connection
close(message_file)  

