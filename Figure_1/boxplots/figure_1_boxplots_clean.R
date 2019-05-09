# Script to make figure 1 boxplot in phylogeny of S. Typhimurium with associated population
# variation in PG (HACS), mean DBS and Invasivness Index for 1st and 3rd population levels

options(stringsAsFactors = FALSE)

library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(gridExtra)

#### LOAD DATA ####

baps.data <- read.table("../baps_data_3lev.txt",
                        sep = "\t", header = TRUE)

HACS.data <- read.table("PG_99.csv",
                          sep = ",", header = TRUE)

invasivness.data <- read.table("mapping_votes.csv",
                         sep = ",", header = TRUE)

dbs.score <- read.table("first_line_core_genome.txt", sep = "\n", header = F, comment.char = "")

#### Clean Data ####

# HACS data
HACS.data$name <- gsub("ALL_REPORTS/report_", "", HACS.data$name)
HACS.data$name <- gsub(".tsv", "", HACS.data$name)
# rownames(HAC.data) <- ariba.data4$name
colnames(HACS.data) <- gsub(".match", "", colnames(HACS.data))
HACS.data[HACS.data == "no"] <- 1
HACS.data[HACS.data == "yes"] <- 0
test.4 <- as.data.frame(lapply(HACS.data, as.numeric))

pg.number <- as.data.frame(cbind(as.data.frame(lapply(HACS.data, as.numeric))[,1],
                                 rowSums(as.data.frame(lapply(HACS.data, as.numeric))
                                         [2:dim(as.data.frame(lapply(HACS.data, as.numeric)))[2]])))
colnames(pg.number) <- c("V2","no.PG")

# DBS data
dbs.split <- sapply(dbs.score, function(x) strsplit(x, ":"))
dbs.mean <- lapply(dbs.split, '[', 3)
dbs.mean <- gsub(" SD", "", dbs.mean)
dbs.mean <- as.numeric(gsub(" ", "", dbs.mean))
dbs.mean <- round(as.numeric(as.character(dbs.mean)),2)
dbs.num <- lapply(dbs.split, '[', 1)
dbs.num <- lapply(dbs.num, function(x) strsplit(x, "_"))
dbs.num <- unlist(dbs.num, recursive = F)
dbs.num <- lapply(dbs.num, '[', 2)
dbs.num <- gsub("/core", "", dbs.num)
dbs.res <- as.data.frame(cbind(dbs.num, dbs.mean))
colnames(dbs.res) <- c("ID","DBS")

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



#### Merge and Prepare Data ####

merge.data <- merge(baps.data[,c(1,4,6)], pg.number, by.x = "SampleID", by.y = "V2")
merge.data <- merge(merge.data, invasivness.data[,c(1,3)], by.x = "SampleID", by.y = "Strain", all = T)
merge.data <- merge(merge.data, dbs.res, by.x = "SampleID", by.y = "ID", all = T)

merge.data <- merge.data[, -1]

colnames(merge.data) <- c("cluster1", "cluster3", "no.PG", "Invasive", "DBS" )

merge.data$cluster3 <- as.character(merge.data$cluster3)
merge.data$cluster1 <- as.character(merge.data$cluster1)
merge.data$DBS <- as.numeric(as.character(merge.data$DBS))

# flatten data for ggplot
flat.data <- melt(merge.data)
# add levels to order data
flat.data$cluster3 <- factor(flat.data$cluster3, levels = c("1","2","3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19"))
flat.data$cluster1 <- factor(flat.data$cluster1, levels = c("1","2"))


#### plot figures ###

p1 <- ggplot(flat.data, aes(x=cluster3,y=value,fill=cluster3)) + geom_boxplot() + facet_wrap(~variable, scales = "free") #+

p2 <- ggplot(flat.data, aes(x=cluster1,y=value,fill=cluster1)) + geom_boxplot() + facet_wrap(~variable, scales = "free") +
    geom_signif(comparisons = list(c("1","2")),
                 map_signif_level=TRUE, test = "t.test")

# arange two plots on grid
plots <- grid.arrange(p2, p1)
#save plot
ggsave("boxplots.pdf",plot = plots,  width = 600, height = 400, units = "mm",
       dpi = 600, bg = "transparent")
