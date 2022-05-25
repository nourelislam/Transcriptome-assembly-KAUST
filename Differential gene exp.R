if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("fission")
BiocManager::install("DESeq2")
BiocManager::install("heatmaps")
BiocManager::install("ComplexHeatmap")
install.packages("pheatmap")
install.packages('dendextend')
install.packages('gplots')
install.packages('pheatmap')
install.packages("dplyr")
install.packages("tidyverse")
install.packages("RColorBrewer")
#########################################################
library(dplyr)
library(tidyverse)
library(tidyverse)
library(readr)
library(ComplexHeatmap)
library(DESeq2)
library(dendextend)
library(ComplexHeatmap)
library("gplots")
library("pheatmap")
library("fission")
library(ggplot2)
library(RColorBrewer)
###########################################################


countData=colnames(simple_transcript_counts)=c("Geneid","33.C.3","33.C.3.2", "33.S.3","33.S.3.2", "33.F.3","33.F.3.2", "33.C.6",
                                          "33.C.6.2", "WT.C.3","WT.C.3.2", "WT.S.3","WT.S.3.2", "33.S.6", "33.S.6.2", "33.F.6",
                                          "33.F.6.2", "WT.F.3","WT.F.3.2", "WT.C.6","WT.C.6.2", "WT.S.6","WT.S.6.2", "WT.F.6","WT.F.6.2")
countData=simple_transcript_counts
rownames(countData) = countData$Geneid

transcripts_only <- read_delim("transcripts_only.txt", ";", escape_double = FALSE, col_names = FALSE, 
                               trim_ws = TRUE)
transcripts_only = transcripts_only[-6]
colnames(transcripts_only) = c("transcripts_id", "gene_id", "ref_transcript_id", "annotation","X")
matching = left_join(countData, transcripts_only, by= c("Geneid" =  "transcripts_id"))

countData=countData[-1]

condition_hours = factor(c("33.cont.3", "33.cont.3", "33.salt.3", "33.salt.3", "33.cold.3", "33.cold.3", "33.cont.6", "33.cont.6", "WT.cont.3", 
                     "WT.cont.3", "WT.salt.3", "WT.salt.3", "33.salt.6","33.salt.6", "33.cold.6", "33.cold.6", "WT.cold.3","WT.cold.3",
                     "WT.cont.6", "WT.cont.6", "WT.salt.6", "WT.salt.6", "WT.cold.6", "WT.cold.6"))

stress = factor(c(rep("unbiotic", 2), rep("abiotic", 2),rep("abiotic", 2), rep("unbiotic", 4), rep("abiotic", 8), rep("unbiotic", 2), rep("abiotic", 4))) 

stress_type = factor(c(rep("33.unbiotic", 2), rep("33.abiotic", 2),rep("33.abiotic", 2), rep("33.unbiotic", 2), rep("WT.unbiotic", 2),
                       rep("WT.abiotic", 2), rep("33.abiotic", 4), rep("WT.abiotic", 2), rep("WT.unbiotic", 2), rep("WT.abiotic", 4))) 
Sample.Type = factor(c(rep("mutant", 8), rep("WT", 4), rep("mutant", 4), rep("WT", 8))) 

env_type = factor(c(rep("33.cont", 2), rep("33.salt", 2), rep("33.cold", 2), rep("33.cont", 2), rep("WT.cont", 2), rep("WT.salt", 2), 
               rep("33.salt", 2),rep("33.cold", 2), rep("WT.cold", 2), rep("WT.cont", 2),rep("WT.salt", 2),rep("WT.cold", 2))) 

env = factor(c(rep("cont", 2), rep("salt", 2), rep("cold", 2), rep("cont", 2), rep("cont", 2), rep("salt", 2), 
               rep("salt", 2),rep("cold", 2), rep("cold", 2), rep("cont", 2),rep("salt", 2),rep("cold", 2)))

X = rownames(colData)
colData = data.frame(Sample.Type, stress, stress_type, env, env_type, condition_hours, X)
rownames(colData) = colnames(countData)
rownames(colData) %in% colnames(countData)

# 1: Comparison of RS33 and WT under control conditions (You can use the 3-hour and 6-hour data under control conditions only).
dds = DESeqDataSetFromMatrix( countData = countData , colData = colData , design = ~ stress_type)
dds.run = DESeq(dds)
### direct results or specifying teh contrast (to make a res object based on two specific conditions/treatment)
res = results(dds.run)
res = results(dds.run, contrast = c("stress_type", "33.unbiotic", "WT.unbiotic"))
# remove nulls
res = res[complete.cases(res), ]
#summary(res.WT.33.cont)
res.df.cont.stress = as.data.frame(res)
plotMA(res, ylim=c(-1,1)) 
summary (res)
res.degs.cont.stress = res.df.cont.stress[res.df.cont.stress$padj< 0.05 & abs(res.df.cont.stress$log2FoldChange)>log2(2),]

res.df.cont.stress$geneid = rownames(res.df.cont.stress)
## sorting and filter the rows ##
filtered_data = arrange(res.df.cont.stress, padj,log2FoldChange)
test = filtered_data[0:500,]
geneids = test$geneid
filter_names = rownames(countData)
## matching ###
#countData$geneid = rownames(countData)
matched = match(geneids, filter_names)
countData$Geneid = rownames(countData)
test = countData$Geneid[matched] ################################## CHANGED
#### indexing ###
top_5hundered= countData[ countData$Geneid %in% test,]
top_5hundered=top_5hundered[-25]
##################################################
######################################
filtered_data = arrange(res.degs.cont.stress, padj,log2FoldChange)
test = filtered_data[0:100,]
test$gene_ids = rownames(test)
Geneids.degs = test$gene_ids
#filter_names = rownames(countData)
## matching ###
#countData$geneid = rownames(countData)
matched.degs = match(Geneids.degs, filter_names)
countData$Geneid = rownames(countData)
idx = countData$Geneid[matched.degs]
#### indexing ###
degs= countData[ countData$Geneid %in% idx,]

countData = countData[-25]

#### transcript annotation ######
degs$transcript.id = rownames(degs)
annotated = left_join(degs, transcripts_only, by= c("transcript.id" =  "transcripts_id"))
annotated[is.na(annotated)] = "u"
rownames(annotated) = degs$transcript.id
annotated = subset(annotated, grepl("Os", annotation))
row_ann = annotated %>% select(Geneid, X) %>% data.frame(row.names = "Geneid")
col_ann = colData %>% select(stress_type) %>% data.frame(row.names = rownames(colData))
degs=degs[-25]
heatmap.degs = pheatmap(degs, scale = "row", annotation_col = col_ann, brewer.pal(5, "YlOrRd"),
                        show_rownames = T, fontsize_row = 5, height = 15,
                        fontsize = 9)

### 33.cont.hours ##
dds = DESeqDataSetFromMatrix( countData = countData , colData = colData , design = ~ condition_hours)
dds.run = DESeq(dds)
### direct results or specifying teh contrast (to make a res object based on two specific conditions/treatment)
res = results(dds.run)
res = results(dds.run, contrast = c("condition_hours", "33.cont.3", "33.cont.6"))
# remove nulls
res = res[complete.cases(res), ]
#summary(res.WT.33.cont)
res.df.cont.33 = as.data.frame(res)
plotMA(res, ylim=c(-1,1)) 
summary (res)
res.degs.cont.33 = res.df.cont.33[res.df.cont.33$padj< 0.05 & abs(res.df.cont.33$log2FoldChange)>log2(2),]

res.df.cont.33$geneid = rownames(res.df.cont.33)
## sorting and filter the rows ##
filtered_data = arrange(res.df.cont.33, padj,log2FoldChange)
test = filtered_data[0:500,]
geneids = test$geneid
filter_names = rownames(countData)
## matching ###
#countData$geneid = rownames(countData)
matched = match(geneids, filter_names)
countData$Geneid = rownames(countData)
test = countData$Geneid[matched] ################################## CHANGED
#### indexing ###
top_5hundered= countData[ countData$Geneid %in% test,]
top_5hundered=top_5hundered[-25]
##################################################
######################################
filtered_data = arrange(res.degs.cont.33, padj,log2FoldChange)
test = filtered_data[0:100,]
test$gene_ids = rownames(test)
Geneids.degs = test$gene_ids
#filter_names = rownames(countData)
## matching ###
#countData$geneid = rownames(countData)
matched.degs = match(Geneids.degs, filter_names)
countData$Geneid = rownames(countData)
idx = countData$Geneid[matched.degs]
#### indexing ###
degs= countData[ countData$Geneid %in% idx,]
col_ann = colData %>% select(condition_hours) %>% data.frame(row.names = rownames(colData))
degs=degs[-25]
heatmap.degs = pheatmap(degs, scale = "row", annotation_col = col_ann, brewer.pal(5, "YlOrRd"),
                        show_rownames = T, fontsize_row = 5, height = 15,
                        fontsize = 9)


#2: Comparison of RS33 and WT under abiotic stress conditions
dds = DESeqDataSetFromMatrix( countData = countData , colData = colData , design = ~stress_type)
dds.run = DESeq(dds)
### direct results or specifying teh contrast (to make a res object based on two specific conditions/treatment)
res=results(dds.run)
res=results(dds.run, contrast = c("stress", "33.abiotic", "WT.abiotic"))
# remove nulls
res=res[complete.cases(res), ]
summary(res)
res.df.stress = as.data.frame(res)
plotMA(res, ylim=c(-1,1)) 
summary (res)
res.degs.stress = res.df.stress[res.df.stress$padj< 0.05 & abs(res.df.stress$log2FoldChange)>log2(2),]

res.df.stress$geneid = rownames(res.df.stress)

degs_annot = degs_annot[-2]
colnames(degs_annot) = "isoform_annot"

#filt = filt[-2]
#rownames(filt) = filt[1]

heatmap_up_500 = pheatmap(top_5hundered, scale = "row")

path= "/media/nour/External Hard Disk/KAUST/splice_site-project/sorted/merged/"
writeLines(rownames(degs), con = "path")

Heatmap(degs, name = "mtcars", 
        row_names_gp = gpar(fontsize = 6.5),
        cluster_rows = color_branches(row_dend, k = 2),
        cluster_columns = color_branches(col_dend, k = 5))
