library(ggplot2)
library(reshape2)


bulk_in <- mymat
bulk_anno <-  s2c
workdir <- "/Users/lab/senguo/projects/sara_rnaseq/deconvolution_v1/"
setwd(workdir)
#
#dir.create(file.path(workdir, "data"))
#dir.create(file.path(workdir, "results"))
###
source("/Users/lab/Downloads/DWLS-master/Deconvolution_functions.R")

##
# load("/Users/lab/Downloads/facs_Skin_seurat_tiss.Robj", verbose=TRUE)
# 
# tiss2 <- UpdateSeuratObject(tiss)

sc_ref <- Read10X("/Users/lab/Downloads/GSE129218/") 

sc_ref2 <- CreateSeuratObject(counts = sc_ref, project = "mice_skin_ref", min.cells = 3, min.features = 200)

###
lables <- read.delim("/Users/lab/Downloads/celltypes_1st_level_9w.txt", sep = "\t",header = F)
lables$V1 <- paste0(lables$V1, "-1")
#
barcodes <- lables[, 1]  # 
cell_types <- lables[, 2]  # A

#
matched_barcodes <-    colnames(sc_ref2) %in% barcodes
#unmatched_barcodes <- barcodes[!matched_barcodes]
sc_ref3 <- sc_ref2[, matched_barcodes]
#
#lables2 <- lables[lables$V1 %in% colnames(sc_ref3), ]
lables2 <- lables[match(colnames(sc_ref3), lables$V1 ), ] ## use *match* to make sure the order is correct, assign of cell types are correct.!!!!!!!!!!!!!!!!!!
#identical(colnames(sc_ref3), lables2$V1) ## double confirm the cell types assginment is correct
sc_ref3$cell.type <- factor(lables2$V2)

#sc_ref3$cell.type <- factor(lables2$V2)
#sc_ref3$cell.type <- factor(cell_types[matched_barcodes])

# load("/Users/lab/Downloads/dataSC.RData", verbose=TRUE)
# load("/Users/lab/Downloads/labels.RData", verbose=TRUE)

# bulk_in_A08 <- bulk_in$A08 ## deconvolute one sample at a time, this is the first sample/column
# names(bulk_in_A08) <- rownames(bulk_in)
# ###
# sc_ref3_mtx <- GetAssayData(sc_ref3)
# sc_ref3_mtx_2 <- as.matrix(sc_ref3_mtx)
# # sc_ref3_mtx_3 <- as.data.frame(sc_ref3_mtx_2)
# label <- as.character(sc_ref3$cell.type)
# #
# Signature<-buildSignatureMatrixMAST(sc_ref3_mtx_2,label,path="results",diff.cutoff=0.5,pval.cutoff=0.01) #### may take 10+ minutes to construct
# ####
# 
# #trim signature and bulk data to contain the same differentially expressed genes
# tr<-trimData(Signature,bulk_in_A08)
# #estimate using dampened weighted least squares
# solDWLS<-solveDampenedWLS(tr$sig,tr$bulk)
# 
# #try running with other estimation methods
# solSVR<-solveSVR(tr$sig,tr$bulk) #nu-SVR
# 
# 
# solOLS<-solveOLS(tr$sig,tr$bulk) #constrained ordinary least squares

###################################################################################
### calculate for all samples of the 3 methods ### this takes several minutes
allCounts_DWLS<-NULL
allCounts_OLS<-NULL
allCounts_SVR<-NULL
for(j in 1:(dim(bulk_in)[2])){
  S<-Signature
  Bulk<-bulk_in[,j]
  names(Bulk)<- rownames(bulk_in) ## modified based on the tutorial , as that input is diff from here !!!!
  Genes<-intersect(rownames(S),names(Bulk))
  B<-Bulk[Genes]
  S<-S[Genes,]
  solOLS<-solveOLS(S,B)
  solDWLS<-solveDampenedWLS(S,B)
  solSVR<-solveSVR(S,B)
  
  allCounts_DWLS<-cbind(allCounts_DWLS,solDWLS)
  allCounts_OLS<-cbind(allCounts_OLS,solOLS)
  allCounts_SVR<-cbind(allCounts_SVR,solSVR)
}
######################################################
colnames(allCounts_DWLS) <- colnames(bulk_in)
allCounts_OLS <- colnames(bulk_in)
allCounts_SVR <- colnames(bulk_in)





# Convert the data frame to long format for plotting
allCounts_DWLS_long <- melt(allCounts_DWLS)

# # Create the stacked bar plot
# ggplot(allCounts_DWLS_long, aes(x = Var2, y = value, fill = Var1)) +
#   geom_bar(stat = "identity") +
#   labs(title = "Stacked Bar Plot of Cell Types by Sample",
#        x = "Sample",
#        y = "Value",
#        fill = "Cell Type") +
#   theme_minimal()
# 
# 
# ###########
# # Define a custom color palette with 8 distinct colors
# custom_colors <- c("#9FE2BF", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#DFFF00")
# 
# # Convert the data frame to long format for plotting
# allCounts_DWLS_long <- melt(allCounts_DWLS)

# 
# pdf("deconvolution_72_mice_rnaseq_v1.pdf",width = 15,height = 9)
# # Create the stacked bar plot with the custom color palette
# ggplot(allCounts_DWLS_long, aes(x = Var2, y = value, fill = Var1)) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = custom_colors) +  # Apply the custom color palette
#   labs(title = "Stacked Bar Plot of Cell Types by Sample",
#        x = "",
#        y = "Ratio",
#        fill = "Cell Type") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Set the font size of x-axis labels
# 
# dev.off()
# write.table(sc_ref3_mtx_2, "sc_mice.txt", quote = F, sep = "\t")
# write.table(lables, "lables_9w.txt", quote = F, sep = "\t", row.names = F)

#####

#save solutions

# save(allCounts_SVR,file="results/allCounts_SVR.RData")
# save(allCounts_DWLS,file="results/allCounts_DWLS.RData")
# save(allCounts_OLS,file="results/allCounts_OLS.RData")


# library(lattice)
# pdf("results/Boxplot_DWLS.pdf",width=10)
# boxplot(allCounts_DWLS[1,1:2],allCounts_DWLS[1,3:5],allCounts_DWLS[1,6:8],allCounts_DWLS[2,1:2],allCounts_DWLS[2,3:5],allCounts_DWLS[2,6:8],allCounts_DWLS[3,1:2],allCounts_DWLS[3,3:5],allCounts_DWLS[3,6:8],colSums(allCounts_DWLS[4:9,1:2]),colSums(allCounts_DWLS[4:9,3:5]),colSums(allCounts_DWLS[4:9,6:8]),names=c("Control","LOF","GOF","Control","LOF","GOF","Control","LOF","GOF","Control","LOF","GOF"),xlab="Condition",ylab="Proportion",ylim=c(0,1))
# stripchart(c(allCounts_DWLS[1,1:2],allCounts_DWLS[1,3:5],allCounts_DWLS[1,6:8],allCounts_DWLS[2,1:2],allCounts_DWLS[2,3:5],allCounts_DWLS[2,6:8],allCounts_DWLS[3,1:2],allCounts_DWLS[3,3:5],allCounts_DWLS[3,6:8],colSums(allCounts_DWLS[4:9,1:2]),colSums(allCounts_DWLS[4:9,3:5]),colSums(allCounts_DWLS[4:9,6:8]))~c(1,1,2,2,2,3,3,3,4,4,5,5,5,6,6,6,7,7,8,8,8,9,9,9,10,10,11,11,11,12,12,12), vertical = TRUE, 
#            method = "jitter", add = TRUE, pch = 20, lwd=2,col = 'blue')
# # text(2,1,"Cycling ISC")
# # text(5,1,"Noncycling ISC")
# # text(8,1,"Transit Amplifying")
# # text(11,1,"Differentiated")
# dev.off()
# 
# pdf("results/Boxplot_OLS.pdf",width=10)
# boxplot(allCounts_OLS[1,1:2],allCounts_OLS[1,3:5],allCounts_OLS[1,6:8],allCounts_OLS[2,1:2],allCounts_OLS[2,3:5],allCounts_OLS[2,6:8],allCounts_OLS[3,1:2],allCounts_OLS[3,3:5],allCounts_OLS[3,6:8],colSums(allCounts_OLS[4:9,1:2]),colSums(allCounts_OLS[4:9,3:5]),colSums(allCounts_OLS[4:9,6:8]),names=c("Control","LOF","GOF","Control","LOF","GOF","Control","LOF","GOF","Control","LOF","GOF"),xlab="Condition",ylab="Proportion",ylim=c(0,1))
# stripchart(c(allCounts_OLS[1,1:2],allCounts_OLS[1,3:5],allCounts_OLS[1,6:8],allCounts_OLS[2,1:2],allCounts_OLS[2,3:5],allCounts_OLS[2,6:8],allCounts_OLS[3,1:2],allCounts_OLS[3,3:5],allCounts_OLS[3,6:8],colSums(allCounts_OLS[4:9,1:2]),colSums(allCounts_OLS[4:9,3:5]),colSums(allCounts_OLS[4:9,6:8]))~c(1,1,2,2,2,3,3,3,4,4,5,5,5,6,6,6,7,7,8,8,8,9,9,9,10,10,11,11,11,12,12,12), vertical = TRUE, 
#            method = "jitter", add = TRUE, pch = 20, lwd=2,col = 'blue')
# # text(2,1,"Cycling ISC")
# # text(5,1,"Noncycling ISC")
# # text(8,1,"Transit Amplifying")
# # text(11,1,"Differentiated")
# dev.off()
# 
# pdf("results/Boxplot_SVR.pdf",width=10)
# boxplot(allCounts_SVR[1,1:2],allCounts_SVR[1,3:5],allCounts_SVR[1,6:8],allCounts_SVR[2,1:2],allCounts_SVR[2,3:5],allCounts_SVR[2,6:8],allCounts_SVR[3,1:2],allCounts_SVR[3,3:5],allCounts_SVR[3,6:8],colSums(allCounts_SVR[4:9,1:2]),colSums(allCounts_SVR[4:9,3:5]),colSums(allCounts_SVR[4:9,6:8]),names=c("Control","LOF","GOF","Control","LOF","GOF","Control","LOF","GOF","Control","LOF","GOF"),xlab="Condition",ylab="Proportion",ylim=c(0,1))
# stripchart(c(allCounts_SVR[1,1:2],allCounts_SVR[1,3:5],allCounts_SVR[1,6:8],allCounts_SVR[2,1:2],allCounts_SVR[2,3:5],allCounts_SVR[2,6:8],allCounts_SVR[3,1:2],allCounts_SVR[3,3:5],allCounts_SVR[3,6:8],colSums(allCounts_SVR[4:9,1:2]),colSums(allCounts_SVR[4:9,3:5]),colSums(allCounts_SVR[4:9,6:8]))~c(1,1,2,2,2,3,3,3,4,4,5,5,5,6,6,6,7,7,8,8,8,9,9,9,10,10,11,11,11,12,12,12), vertical = TRUE, 
#            method = "jitter", add = TRUE, pch = 20, lwd=2,col = 'blue')
# # text(2,1,"Cycling ISC")
# # text(5,1,"Noncycling ISC")
# # text(8,1,"Transit Amplifying")
# # text(11,1,"Differentiated")
# dev.off()


##############################################################
# allCounts_DWLS2 <- as.data.frame(allCounts_DWLS)
# 
# EPI_ANA_sums <- colSums(allCounts_DWLS2[c("EPI", "ANA"), ])
# 
# bulk_in_kera <- bulk_in
# 
# for (i in 1:ncol(bulk_in)) {
#   bulk_in_kera[, i] <- bulk_in[, i] * EPI_ANA_sums[i]
# }
# 
# head(bulk_in_kera)
# ##############################################################
# 
# 
# WriteXLS(c("bulk_in_kera"),"",row.names=F,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)


###############################################################
bulk_anno2 <- bulk_anno %>% as.data.frame() %>% dplyr::select(-path)
bulk_anno2$Sample <- rownames(bulk_anno2)

####
# pdf("deconvolution_72_mice_rnaseq_v2.pdf",width = 15,height = 9)
# # Create the stacked bar plot with the custom color palette
# ggplot(allCounts_DWLS_long, aes(x = Var2, y = value, fill = Var1)) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = custom_colors) +  # Apply the custom color palette
#   labs(title = "Stacked Bar Plot of Cell Types by Sample",
#        x = "",
#        y = "Ratio",
#        fill = "Cell Type") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Set the font size of x-axis labels
# 
# dev.off()

###########################################################################
##########################################################################
################################################################################
# # Transform the allCounts_DWLS data frame to a long format
# counts_long <- allCounts_DWLS %>%
#   as.data.frame() %>%
#   rownames_to_column("CellType") %>%
#   gather(key = "Sample", value = "Count", -CellType)
# 
# # Merge with bulk_anno2 to get the annotation details
# merged_data <- merge(counts_long, bulk_anno2, by.x = "Sample", by.y = "Sample")
# 
# # Set factor levels for ordering
# merged_data$Time <- factor(merged_data$Time, levels = c("t0", "t24", "t30"))
# merged_data$Group <- factor(merged_data$Group, levels = c("Ctrl_0", "Ctrl_24", "Ctrl_30", "Blue", "UVB", "UVC"))
# merged_data$Genotype <- factor(merged_data$Genotype, levels = c("BL6", "NZM"))
# 
# # Plot
# ggplot(merged_data, aes(x = Sample, y = Count, fill = CellType)) +
#   geom_bar(stat = "identity") +
#   facet_wrap(~Time, scales = "free_x") +
#   labs(title = "Stacked Bar Plot of Cell Type Counts",
#        x = "Sample",
#        y = "Count",
#        fill = "Cell Type") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   geom_vline(data = merged_data %>% group_by(Time, Group) %>% summarize(n = n()),
#              aes(xintercept = cumsum(n) - 0.5),
#              linetype = "dashed")

###################################################################################################
#####
#######################################################################################
# # Transform the allCounts_DWLS data frame to a long format
# counts_long <- allCounts_DWLS %>%
#   as.data.frame() %>%
#   rownames_to_column("CellType") %>%
#   gather(key = "Sample", value = "Count", -CellType)
# 
# # Merge with bulk_anno2 to get the annotation details
# merged_data <- merge(counts_long, bulk_anno2, by.x = "Sample", by.y = "Sample")
# 
# # Set factor levels for ordering
# merged_data$Time <- factor(merged_data$Time, levels = c("t0", "t24", "t30"))
# merged_data$Group <- factor(merged_data$Group, levels = c("Ctrl_0", "Ctrl_24", "Ctrl_30", "Blue", "UVB", "UVC"))
# merged_data$Genotype <- factor(merged_data$Genotype, levels = c("BL6", "NZM"))
# 
# # Create a unique ordering based on Time, Group, and Genotype
# unique_order <- merged_data %>%
#   arrange(Time, Group, Genotype) %>%
#   distinct(Sample, .keep_all = TRUE) %>%
#   pull(Sample)
# 
# # Assign the unique order to the Sample factor
# merged_data <- merged_data %>%
#   mutate(Sample = factor(Sample, levels = unique_order))
# 
# # Calculate the number of samples in each time group for adjusting bar width
# time_group_counts <- merged_data %>% group_by(Time) %>% summarize(n = n_distinct(Sample), .groups = 'drop')
# 
# # Plot
# ggplot(merged_data, aes(x = Sample, y = Count, fill = CellType)) +
#   geom_bar(stat = "identity", width = 0.7) +
#   facet_wrap(~Time, scales = "free_x") +
#   labs(title = "Stacked Bar Plot of Cell Type Counts",
#        x = "",
#        y = "Count",
#        fill = "Cell Type") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   geom_vline(data = merged_data %>% group_by(Time, Group) %>% summarize(n = n(), .groups = 'drop'),
#              aes(xintercept = cumsum(n) - 0.5),
#              linetype = "dashed") +
#   scale_x_discrete(labels = function(x) {
#     ifelse(merged_data$Time[which(levels(merged_data$Sample) == x)[1]] == "t30", x, "")
#   }) +
#   theme(strip.background = element_blank())
# Load necessary libraries
# library(ggplot2)
# library(dplyr)
# library(tidyr)
# library(ggforce)
# 
# # Assuming allCounts_DWLS and bulk_anno2 data frames are already loaded
# 
# # Transpose allCounts_DWLS and convert it to a long format
# allCounts_DWLS_long <- allCounts_DWLS %>%
#   as.data.frame() %>%
#   rownames_to_column(var = "CellType") %>%
#   gather(key = "Sample", value = "Count", -CellType)
# 
# # Join with bulk_anno2 to get additional annotations
# allCounts_DWLS_long <- allCounts_DWLS_long %>%
#   inner_join(bulk_anno2, by = c("Sample" = "Name"))
# 
# # Sort by Time and Genotype
# allCounts_DWLS_long <- allCounts_DWLS_long %>%
#   arrange(Time, Genotype)
# 
# # Assign relative widths to each time point
# time_levels <- unique(allCounts_DWLS_long$Time)
# widths <- ifelse(time_levels == "t30", 3, 1) # Adjust width for t30
# 
# allCounts_DWLS_long$Time <- factor(allCounts_DWLS_long$Time, 
#                                    levels = c("t0", "t24", "t30"),
#                                    labels = c("T0", "T24", "T30"))
# 
# # Plot the stacked bar plot with adjusted facet widths
# ggplot(allCounts_DWLS_long, aes(x = reorder(Sample, Time), y = Count, fill = CellType)) +
#   geom_bar(stat = "identity") +
#   ggforce::facet_col(~ Time, scales = "free_x", space = "free") +
#   theme_minimal() +
#   labs(title = "Stacked Bar Plot of Cell Type Counts",
#        x = "",
#        y = "Predicted Ratio",
#        fill = "Cell Type") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size =12),
#         strip.text = element_text(size = 12),  # Adjusted font size for facet strip labels
#         axis.title.y = element_text(size = 12),
#         axis.text.y = element_text(size = 12))



# # Load necessary libraries
# library(ggplot2)
# library(dplyr)
# library(tidyr)
# library(ggforce)
# 
# # Assuming allCounts_DWLS and bulk_anno2 data frames are already loaded
# 
# # Transpose allCounts_DWLS and convert it to a long format
# allCounts_DWLS_long <- allCounts_DWLS %>%
#   as.data.frame() %>%
#   rownames_to_column(var = "CellType") %>%
#   gather(key = "Sample", value = "Count", -CellType)
# 
# # Join with bulk_anno2 to get additional annotations
# allCounts_DWLS_long <- allCounts_DWLS_long %>%
#   inner_join(bulk_anno2, by = c("Sample" = "Name"))
# 
# # Create a new sorting column
# allCounts_DWLS_long <- allCounts_DWLS_long %>%
#   mutate(SortOrder = case_when(
#     Time %in% c("t0", "t24") ~ paste(Genotype, Sample, sep = "_"),
#     Time == "t30" ~ paste(Genotype, Group, Sample, sep = "_")
#   ))
# 
# # Convert SortOrder to a factor to ensure the correct order
# allCounts_DWLS_long$SortOrder <- factor(allCounts_DWLS_long$SortOrder, 
#                                         levels = unique(allCounts_DWLS_long$SortOrder[order(allCounts_DWLS_long$Time, allCounts_DWLS_long$Genotype, allCounts_DWLS_long$Group, allCounts_DWLS_long$Sample)]))
# 
# # Rename the levels of Time factor
# allCounts_DWLS_long$Time <- factor(allCounts_DWLS_long$Time, 
#                                    levels = c("t0", "t24", "t30"),
#                                    labels = c("T0", "T24", "T30"))
# 
# # Plot the stacked bar plot with adjusted facet widths
# ggplot(allCounts_DWLS_long, aes(x = SortOrder, y = Count, fill = CellType)) +
#   geom_bar(stat = "identity") +
#   ggforce::facet_col(~ Time, scales = "free_x", space = "free") +
#   theme_minimal() +
#   labs(title = "Stacked Bar Plot of Cell Type Counts",
#        x = "",
#        y = "Predicted Ratio",
#        fill = "Cell Type") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size =12),
#         strip.text = element_text(size = 12),  # Adjusted font size for facet strip labels
#         axis.title.y = element_text(size = 12),
#         axis.text.y = element_text(size = 12))
##############################################################################################################
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggforce)

# Assuming allCounts_DWLS and bulk_anno2 data frames are already loaded

# Transpose allCounts_DWLS and convert it to a long format
allCounts_DWLS_long <- allCounts_DWLS %>%
  as.data.frame() %>%
  rownames_to_column(var = "CellType") %>%
  gather(key = "Sample", value = "Count", -CellType)

# Join with bulk_anno2 to get additional annotations
allCounts_DWLS_long <- allCounts_DWLS_long %>%
  inner_join(bulk_anno2, by = c("Sample" = "Name"))

# Define the order of groups within t30
group_order <- c("Ctrl_30", "Blue", "UVB", "UVC")

# Create a new sorting column with specific order for t30 groups
allCounts_DWLS_long <- allCounts_DWLS_long %>%
  mutate(Group = factor(Group, levels = group_order),
         SortOrder = case_when(
           Time %in% c("t0", "t24") ~ paste(Genotype, Sample, sep = "_"),
           Time == "t30" ~ paste(Genotype, Group, Sample, sep = "_")
         ))

# Convert SortOrder to a factor to ensure the correct order
allCounts_DWLS_long$SortOrder <- factor(allCounts_DWLS_long$SortOrder, 
                                        levels = unique(allCounts_DWLS_long$SortOrder[order(allCounts_DWLS_long$Time, allCounts_DWLS_long$Genotype, allCounts_DWLS_long$Group, allCounts_DWLS_long$Sample)]))

# Rename the levels of Time factor
allCounts_DWLS_long$Time <- factor(allCounts_DWLS_long$Time, 
                                   levels = c("t0", "t24", "t30"),
                                   labels = c("T0", "T24", "T30"))


custom_colors2 <- c("#A94F1A", "#4E8C65", "#6359A2", "#B52177", "#ff7f0e", "#C69E24",  "#666666", "#6A9623")


pdf("sarah_rna_UV_72_deconvolution_v2_2_correct_sample_info.pdf",width = 16,height = 9)
# Plot the stacked bar plot with adjusted facet widths
ggplot(allCounts_DWLS_long, aes(x = SortOrder, y = Count, fill = CellType)) +
  geom_bar(stat = "identity") + scale_fill_manual(values = custom_colors2) + 
  ggforce::facet_col(~ Time, scales = "free_x", space = "free") +
  theme_minimal() +
  labs(title = "Stacked Bar Plot of Cell Type Counts",
       x = "",
       y = "Predicted Ratio",
       fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        strip.text = element_text(size = 12),  # Adjusted font size for facet strip labels
        axis.title.y = element_text(size = 13),
        axis.text.y = element_text(size = 12),
        plot.margin = unit(c(0.8, 0.8, 0.8, 1.1), "cm"))
dev.off()
################################################
####
########################




