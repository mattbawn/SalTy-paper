# Script to make figure 1 in phylogeny of S. Typhimurium with associated population
# structure and presence / absence of Plasmid, AMR and HAC genes

options(stringsAsFactors = FALSE)

library(ggtree)
library(phangorn)
library(RColorBrewer)
library(ggplot2)

### LOAD DATA ###

# tree
my.tree <- read.tree("Typhimurium-tree.tre")

# Plasmid, AMR and VF gene data from ARIBA
ariba.data <- read.table("ariba-data.tsv",
                         sep = "\t", header = TRUE)

# HACS data
srst2.data <- read.table("Supplementary_table_4.tsv",
                         sep = "\t", header = TRUE)
# switch HAC and WT in lpfD
srst2.data$lpfD[srst2.data$lpfD == "HAC"] <- "x"
srst2.data$lpfD[srst2.data$lpfD == "WT"] <- "HAC"
srst2.data$lpfD[srst2.data$lpfD == "x"] <- "WT"

# Sample metadata
sample.data <- read.csv("samples-overview-representative_host.csv", stringsAsFactors = FALSE, sep = ",")

# BAPS data
baps.data <- read.table("baps_data_3lev.txt",
                        sep = "\t", header = TRUE)

### Clean data ###

# clean sample data
sample.data$source[sample.data$source == "turkey" | 
                     sample.data$source == "chicken" ] <- "poultry"
sample.data$source[sample.data$source == "duck" | 
                     sample.data$source == "duck " ] <- "duck"
sample.data$source[sample.data$source == "cattle" | 
                     sample.data$source == "cattle " ] <- "cattle"
sample.data$source[sample.data$source == "feed" | 
                     sample.data$source == "environment" ] <- "environment"
sample.data$source[sample.data$source == "dog" | sample.data$source == "cat" |  
                     sample.data$source == "sheep" |  sample.data$source == "horse" ] <- "other"
sample.data$source[sample.data$source == "pigeon" | 
                     sample.data$source == "bird" ] <- "wild bird"
sample.data$source[sample.data$source == "pig" | sample.data$source == "pig " ] <- "pig"

# clean ariba data
ariba.data$name <- as.character(ariba.data$name)
colnames(ariba.data) <- gsub(".match", "", colnames(ariba.data))
ariba.data[ariba.data == "Absent"] <- ""

# Rownames have to correspond to tip labels
rownames(ariba.data) <- ariba.data[,1]
ariba.data <- ariba.data[,-1]

# Format Tree
# Remove Reference from tree
my.tree <- drop.tip(my.tree, "Reference")
# reroot tree based on Heidelberg as outgroup
my.tree <- reroot(my.tree, 141)

# Format Tree tip labels
my.tree$tip.label[grep("mysnps_rep_", my.tree$tip.label)] <- 
  sapply(my.tree$tip.label[grep("mysnps_rep_", my.tree$tip.label)], 
         function(x) as.numeric(gsub("mysnps_rep_","", x))+120)

my.tree$tip.label <- sapply(my.tree$tip.label, function(x) as.numeric(gsub("mysnps_","", x)))
my.tree$tip.label <- as.character(my.tree$tip.label)

# BAPS data
colnames(baps.data) <- c("SampleID","Row", "X1.layer.clusterID", "level1","level2","level3" )

# Change sample names 
baps.data$SampleID[grep("mysnps_rep_", baps.data$SampleID)] <- 
  sapply(baps.data$SampleID[grep("mysnps_rep_", baps.data$SampleID)], 
         function(x) as.numeric(gsub("mysnps_rep_","", x))+120)

# Change sample names 
baps.data$SampleID[grep("mysnps_", baps.data$SampleID)] <- 
  sapply(baps.data$SampleID[grep("mysnps_", baps.data$SampleID)], 
         function(x) as.numeric(gsub("mysnps_","", x)))


### Merge Input Data ####

merge.data <- as.data.frame(merge(sample.data, baps.data, by.x = "V2", by.y = "SampleID", all = TRUE))


# Rownames have to correspond to tip labels
rownames(merge.data) <- merge.data$V2
merge.data$number <- merge.data$V2
merge.data$level3 <- as.character(merge.data$level3)

# short ariba and HACS
short.data <- as.data.frame(merge(ariba.data[,1:44], srst2.data[,-1],by="row.names",all=TRUE))
rownames(short.data) <- short.data[,1]
short.data <- short.data[,-c(1,dim(short.data)[2])]


# colour data
short.data[short.data == ""] <- NA
cols.1=unlist(short.data)

cols.1[cols.1 == "WT"] <- "lightgrey"
cols.1[cols.1 == "WT?"] <- "grey"
cols.1[cols.1 == "HAC"] <- "deeppink"
cols.1[cols.1 == "HAC?"] <- "deeppink2"
cols.1[cols.1 == "HAC2"] <- "orange"
cols.1[cols.1 == "HAC3"] <- "darkorange"
cols.1[cols.1 == "Resfinder"] <- "green3"
cols.1[cols.1 == "Plasmidfinder"] <- "red3"

names(cols.1) <- unlist(short.data)

### plot ###

par(bg="transparent")
p <- ggtree(my.tree, size=1)
p <- p + geom_tiplab(size=2, align = TRUE)

# add gene data to tree
p1 <- gheatmap(p, short.data,  offset = 0.02, width=3, colnames = F) %>%
  scale_x_ggtree() + scale_fill_manual(values=cols.1)

# add metadata to tree and color tip labels
p2 <- p1 %<+% merge.data + geom_tiplab(aes(label=type, color=level3), size=5, align =TRUE , show.legend = FALSE)


df <- get_heatmap_column_position(p2, by="top")
p3 <- p2 + geom_text(data = df, aes(x,y+2, label=label), size=1.5, angle=90 ) + geom_tippoint(aes(color=source),size=3)

# # ggsave("rerroted-tree.pdf",plot = p5,  width = 800, height = 500, units = "mm",
# #         bg = "transparent")
# 
