### Loading libraries ####

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(ggpubr)


# ### Preparing the tables ### 

# loading RNAseq information of treatments

zebra_general <- read.table('data/Conectivity_Zebrafish.txt', sep = "\t", stringsAsFactors = F, header = T)

amphi_general <- read.table('data/Conectivity_amphi.txt', sep = "\t", stringsAsFactors = F, header = T)

#Loading the connectivity of genes that respond to ATAC and RNAseq

zebra_atac_rna <- read.table('data/Conectivity_Zebra_ATAC_RNA_restrictive.txt', sep = "\t", stringsAsFactors = F, header = T)

amphi_atac_rna <- read.table('data/Conectivity_Amphi_ATAC_vs_RNASeq.txt', sep = "\t", stringsAsFactors = F, header = T)

# Mending the colnames, so every table has the same header

colnames(amphi_general) <- colnames(zebra_general)

colnames(zebra_atac_rna) <- colnames(zebra_general)

colnames(amphi_atac_rna) <- colnames(zebra_general)

merge <- rbind(zebra_general, amphi_general, zebra_atac_rna, amphi_atac_rna)

## Generation of a table that has all the information that will be needed by ggplot2

merge$Cols <- factor(rep(c('Zebrafish', 'Amphioxus', 'Zebrafish', 'Amphioxus'),
                            times = c(nrow(zebra_general), nrow(amphi_general), nrow(zebra_atac_rna), nrow(amphi_atac_rna))))

merge$Rows <- factor(rep(c('RNASeq', 'RNASeq', 'ATAC & RNASeq', 'ATAC & RNASeq'),
                         times = c(nrow(zebra_general), nrow(amphi_general), nrow(zebra_atac_rna), nrow(amphi_atac_rna))))

merge$Species <- factor(rep(c('Zebrafish RNASeq', 'Amphioxus RNASeq', 'Zebrafish ATAC & RNASeq', 'Amphioxus ATAC & RNASeq'),
                            times = c(nrow(zebra_general), nrow(amphi_general), nrow(zebra_atac_rna), nrow(amphi_atac_rna))))

colnames(merge)[6] <- "Connectivity"

write.table(merge, file = 'outputs/Amphi_zebra_Connectivity.txt', sep = '\t', quote = F, row.names = F)

## Connectivity Zebrafish & Amphi General####

data <- read.table('outputs/Amphi_zebra_Connectivity.txt', sep = "\t", stringsAsFactors = F, header = T)

data$Rows <- relevel(x = factor(data$Rows), ref = 'RNASeq')

myColors <- brewer.pal(4, "Paired")
names(myColors) <- c(1,2,3,4)

# making the plot of connectivity using the previous data.

ggplot(data, aes(x= factor(Connectivity),group=factor(Species))) + # in the x we have connectivity and the two species
  geom_bar(aes(y= ..prop.., fill = factor(Cols)), stat="count", position = 'dodge2')+ 
  facet_wrap(~Rows, ncol =  2, scales='free') + # we have the info of only RNAseq and also RNAseq + ATACseq
  theme_classic()  +
  scale_y_continuous(limits=c(0,1),labels = scales::percent) +
  theme(plot.title = element_text(hjust = 0.5, size = 12), plot.subtitle = element_text(hjust = 0.5), 
        axis.text = element_text(size = 12), axis.title=element_text(size=12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 12),
        strip.text.x = element_text(size = 12), panel.spacing = unit(8, "mm"), strip.text.y = element_text(size = 12)) +
  ggtitle('Connectivity of Responsive Genes') + 
  guides(fill = guide_legend(title = 'Species')) + xlab('Connectivity') + ylab('Relative Frequences')+
  geom_text(aes(label=scales::percent(..prop..), 
                y=..prop..), stat="count", vjust=-.5,position = position_dodge(width = 1)) + 
  scale_fill_manual(values = c('#1b9e77', '#377eb8'))
  

ggsave('outputs/Connectivity_amphi_zebra.pdf', width = 210 ,units = 'mm')

ggsave('outputs/Connectivity_amphi_zebra.png', width = 210 ,units = 'mm')


#### Statistics of proportions ####

# In Rnaseq

data_rnaseq <- data %>% filter(Rows == "RNASeq")

for (connectivity in 1:4) {
  
  print(paste("Connectivity:", connectivity))
  
  amphi_genes <- nrow(data_rnaseq %>% filter(Connectivity == connectivity, Cols == "Amphioxus"))
  amphi_total <- nrow(data_rnaseq %>% filter(Cols == "Amphioxus"))
  
  zebra_genes <- nrow(data_rnaseq %>% filter(Connectivity == connectivity, Cols == "Zebrafish"))
  zebra_total <- nrow(data_rnaseq %>% filter(Cols == "Zebrafish"))
  
  print(prop.test(x = c(amphi_genes, zebra_genes), n = c(amphi_total, zebra_total)))
  print(prop.test(x = c(amphi_genes, zebra_genes), n = c(amphi_total, zebra_total))$p.value)
}

#In ATAC & RNaseq

data_ATAC_RNA <- data %>% filter(Rows == "ATAC & RNASeq")

for (connectivity in 1:4) {
  
  print(paste("Connectivity:", connectivity))
  
  amphi_genes <- nrow(data_ATAC_RNA %>% filter(Connectivity == connectivity, Cols == "Amphioxus"))
  amphi_total <- nrow(data_ATAC_RNA %>% filter(Cols == "Amphioxus"))
  
  zebra_genes <- nrow(data_ATAC_RNA %>% filter(Connectivity == connectivity, Cols == "Zebrafish"))
  zebra_total <- nrow(data_ATAC_RNA %>% filter(Cols == "Zebrafish"))
  
  print(prop.test(x = c(amphi_genes, zebra_genes), n = c(amphi_total, zebra_total)))
}



######### Using Number of copys ########

# ## load the data
 ohnologs <- read.delim('data/ohnolog_list.tsv', stringsAsFactors = F)

 # Generate a intermediate variable with all Bla gene names

 ohnologs_bla <- ohnologs$Bla

 names(ohnologs_bla) <- ohnologs$Dre

 ohnologs_bla <- unique(ohnologs_bla)

 # this for loop will generate 2 vectors:
 #     -bla_and_dre: Is a char vector containing all identifiers, ordered by ohnolog fam
 #     -ohn_fam: a char vector of the same length than bla_and_dre, which will save the info of to which ohnolog family the gene belongs to.

 bla_and_dre <- NULL
 ohn_fam <- NULL

 for (i in ohnologs_bla) {

   bla_and_dre_i <- c(i, ohnologs[ohnologs$Bla == i, "Dre"])

   bla_and_dre <- c(bla_and_dre, bla_and_dre_i)

   ohn_fam_i <- paste0('ohn_fam_',i)

   ohn_fam <- c(ohn_fam,rep(ohn_fam_i, length(bla_and_dre_i)))

 }

 # we check that we dont have any duplicates in our set of ID
 any(duplicated(bla_and_dre))

 # We take the info (mmcat and drecat) about Bla and Dre genes and fix the data frames

 ohnologs_Bla_only <- ohnologs %>%
   select(Bla,mm_cat, dre_cat) %>%
   filter(!duplicated(Bla))

 colnames(ohnologs_Bla_only) <- c("tss_id", "mm_cat",  "dre_cat")
 rownames(ohnologs_Bla_only) <-  as.character(ohnologs_Bla_only$tss_id)

 ohnologs_dre_only <- ohnologs %>%
   select(Dre, mm_cat, dre_cat) %>%
   filter(!duplicated(Dre))

 colnames(ohnologs_dre_only) <- c("tss_id", "mm_cat",  "dre_cat")
 rownames(ohnologs_dre_only) <- as.character(ohnologs_dre_only$tss_id)

 ## We define what will be the base data.frame. We will use this data structure for subsequent analysis

 ohnologs_fixed <- data.frame(tss_id = as.character(bla_and_dre), ohn_fam = ohn_fam)

 rownames(ohnologs_fixed) <- as.character(ohnologs_fixed$tss_id)

 ## We merge the base data.frame with the previous ones, so now we have a structure that can be worked and also the onolog info

 ohnologs_fixed <- full_join(ohnologs_fixed, ohnologs_dre_only, "tss_id") %>%
   full_join(ohnologs_Bla_only, "tss_id")

 ohnologs_fixed$mm_cat.x[is.na(ohnologs_fixed$mm_cat.x)] <- ohnologs_fixed[!is.na(ohnologs_fixed$mm_cat.y),"mm_cat.y"]

 ohnologs_fixed$dre_cat.x[is.na(ohnologs_fixed$dre_cat.x)] <- ohnologs_fixed[!is.na(ohnologs_fixed$dre_cat.y),"dre_cat.y"]

 ohnologs_fixed <- ohnologs_fixed[,c(1:4)]

 colnames(ohnologs_fixed) <- c("tss_id","ohn_fam", "mm_cat",  "dre_cat")

 ## Just a reconversion of the data frame to recycle old code

 boundaries <- ohnologs_fixed

 rownames(boundaries) <- boundaries$tss_id

 ## Loading Connectivity data into the copy data ##

 ## loading conectivity for each species and select only the genes and general conectivity

 zebra_general <- read.table('data/Conectivity_Zebrafish.txt', sep = "\t", stringsAsFactors = F, header = T)

 amphi_general <- read.table('data/Conectivity_amphi.txt', sep = "\t", stringsAsFactors = F, header = T)


 colnames(zebra_general) <- c("tss_id", "BIO",      "RA",         "SB505",      "SU5402" ,    "Conectivity")
 colnames(amphi_general) <- c("tss_id", "BIO",      "RA",         "SB505",      "SU5402" ,    "Conectivity")

 amphi_general <- amphi_general %>% select(tss_id, Conectivity)

 zebra_general <- zebra_general %>% select(tss_id, Conectivity)

 ## Merge of the copy info(boundaries) to connectivity info. We do the merging regarding of the tss_id column
 # We need to fix also the connectivity column, in order to have just one

 boundaries_conect <- boundaries %>% full_join(zebra_general, by = "tss_id") %>% left_join(amphi_general, by = "tss_id") %>% filter(!is.na(dre_cat))

 boundaries_conect$Conectivity.x[is.na(boundaries_conect$Conectivity.x)] <- 0

 boundaries_conect$Conectivity.y[is.na(boundaries_conect$Conectivity.y)] <- 0

 ## Finally we take only those columns of interest

 boundaries_conect <- boundaries_conect %>%
   mutate(Connectivity = Conectivity.x + Conectivity.y) %>%
   select(tss_id,ohn_fam, mm_cat, dre_cat, Connectivity)

 write.table(boundaries_conect, file = 'outputs/ohnologs_connected.tsv', sep = '\t', col.names = T, row.names = F, quote = F)

 
 # Now lets plot this data and see the effect of copy retention in connectivity

boundaries_filtered <-  read.table('outputs/ohnologs_connected.tsv', sep = "\t", stringsAsFactors = F, header = T)

# genes of category 1:4 or more will be considered the same

boundaries_filtered[boundaries_filtered$dre_cat > 4, "dre_cat"] <- 4


boundaries_filtered %>% 
ggplot(aes(x=as.factor(dre_cat))) +
  geom_bar(position = "fill", aes(y= ..count../sum(..count..), fill =factor(Connectivity, levels = c("4","3","2","1","0"))), width = 0.7) + 
  scale_fill_brewer(palette = 'Dark2', direction = -1) +
  scale_y_continuous(labels = scales::percent_format()) + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 12), plot.subtitle = element_text(hjust = 0.5), 
        axis.text = element_text(size = 12), axis.title=element_text(size=12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 12),
        strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12)) + 
  scale_x_discrete(labels = c("1:1", "1:2", "1:3", "1:4 or more")) +
  xlab('Gene category') + ylab('Relative Frequences') + guides(fill = guide_legend(title = 'Connectivity'))

barplot_connect_by_gene_category <- boundaries_filtered %>% 
  ggplot(aes(x=as.factor(dre_cat))) +
  geom_bar(position = "fill", aes(y= ..count../sum(..count..), fill =factor(Connectivity, levels = c("4","3","2","1","0"))), width = 0.7) + 
  scale_fill_brewer(palette = 'Dark2', direction = -1) +
  scale_y_continuous(labels = scales::percent_format()) + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 12), plot.subtitle = element_text(hjust = 0.5), 
        axis.text = element_text(size = 12), axis.title=element_text(size=12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 12),
        strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12)) + 
  scale_x_discrete(labels = c("1:1", "1:2", "1:3", "1:4 or more")) +
  xlab('Gene category') + ylab('Relative Frequences') + guides(fill = guide_legend(title = 'Connectivity'))


ggsave('outputs/Connect_by_copy_danrer.png', width = 210 ,units = 'mm')

ggsave('outputs/Connect_by_copy_danrer.pdf', width = 210 ,units = 'mm')
