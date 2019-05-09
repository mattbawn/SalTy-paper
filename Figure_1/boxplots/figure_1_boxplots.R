
options(stringsAsFactors = FALSE)


ariba.data <- read.table("resfinder.csv",
                         sep = ",", header = TRUE)
ariba.data2 <- read.table("VFDB_full_DT2s.csv",
                         sep = ",", header = TRUE)

ariba.data3 <- read.table("DBS.csv",
                          sep = ",", header = TRUE)

ariba.data4 <- read.table("PG_99.csv",
                          sep = ",", header = TRUE)
ariba.data$name <- gsub("ALL_REPORTS/report_", "", ariba.data$name)
ariba.data$name <- gsub(".tsv", "", ariba.data$name)
ariba.data$name <- gsub("output_resfinder/", "", ariba.data$name)

ariba.data2$name <- gsub("ALL_REPORTS/report_", "", ariba.data2$name)
ariba.data2$name <- gsub(".tsv", "", ariba.data2$name)
ariba.data2$name <- gsub("output_vfdbfull/", "", ariba.data2$name)

ariba.data3$name <- gsub("ALL_REPORTS/report_", "", ariba.data3$name)
ariba.data3$name <- gsub(".tsv", "", ariba.data3$name)

ariba.data4$name <- gsub("ALL_REPORTS/report_", "", ariba.data4$name)
ariba.data4$name <- gsub(".tsv", "", ariba.data4$name)


snippy.snps <- read.table("SNIPPY_SNPS.vcf", header = TRUE, sep = "\t", comment.char = "", skip = 3)
snippy.snps <- snippy.snps[,-(1:10)]
no.snps <- colSums(snippy.snps)
snps.lab <- names(no.snps)

srst2.data <- read.table("mapping_votes.csv",
                         sep = ",", header = TRUE)
# srst2.data$Strain <- lapply(strsplit(srst2.data$Strain, "_"), '[',1)

srst2.data <- as.data.frame(srst2.data)

# srst2.data$Strain <- unlist(srst2.data$Strain)

reps <- grep("mysnps_rep_", snps.lab)

snps.lab[reps] <- sapply(snps.lab[reps], function(x) as.numeric(gsub("mysnps_rep_","", x))+120)

snps.lab <- sapply(snps.lab, function(x) as.numeric(gsub("mysnps_","", x)))

snps <- as.data.frame(cbind(snps.lab, no.snps))
colnames(snps) <- c("ID", "snps")

rownames(ariba.data) <- ariba.data$name
rownames(ariba.data2) <- ariba.data2$name
rownames(ariba.data3) <- ariba.data3$name
rownames(ariba.data4) <- ariba.data4$name


colnames(ariba.data) <- gsub(".match", "", colnames(ariba.data))
colnames(ariba.data2) <- gsub(".match", "", colnames(ariba.data2))
colnames(ariba.data3) <- gsub(".match", "", colnames(ariba.data3))
colnames(ariba.data4) <- gsub(".match", "", colnames(ariba.data4))

sample.data <- read.csv("../samples-overview-representative_host.csv", stringsAsFactors = FALSE, sep = ",")

gifsy1 <- read.table("Gifsy1.txt", sep = "\t", header = T)
gifsy2 <- read.table("Gifsy2.txt", sep = "\t", header = T)

dbs.score <- read.table("first_line_core_genome.txt", sep = "\n", header = F, comment.char = "")
dbs.split <- sapply(dbs.score, function(x) strsplit(x, ":"))
dbs.mean <- lapply(dbs.split, '[', 3)
dbs.mean <- gsub(" SD", "", dbs.mean)
dbs.mean <- as.numeric(gsub(" ", "", dbs.mean))
dbs.num <- lapply(dbs.split, '[', 1)
dbs.num <- lapply(dbs.num, function(x) strsplit(x, "_"))
dbs.num <- unlist(dbs.num, recursive = F)
dbs.num <- lapply(dbs.num, '[', 2)
# dbs.num <- gsub("/SL1344", "", dbs.num)
dbs.num <- gsub("/core", "", dbs.num)
dbs.res <- as.data.frame(cbind(dbs.num, dbs.mean))

colnames(dbs.res) <- c("ID","DBS")

gifsy1$ID <- gsub("mysnps_Gifsy1_", "",gifsy1$ID)
gifsy2$ID <- gsub("mysnps_Gifsy2_", "",gifsy2$ID)

gifsy1 <- gifsy1[,-c(2,3)]
gifsy2 <- gifsy2[,-c(2,3)]

colnames(gifsy1) <- c("ID", "Gifsy1")
colnames(gifsy2) <- c("ID", "Gifsy2")

Gifsys <- merge(gifsy1, gifsy2, by.x = "ID", by.y = "ID")
Gifsys_DBS <- merge(Gifsys, dbs.res, by.x = "ID", by.y = "ID")
Gifsys_DBS_snps <- merge(Gifsys_DBS, snps, by.x = "ID", by.y = "ID")
# Remove 24, 48 and 49 from data

short.sample.data <- sample.data[-c(24,48,49),]



ariba.data[ariba.data == "no"] <- 0
ariba.data[ariba.data == "yes"] <- 1
ariba.data2[ariba.data2 == "no"] <- 0
ariba.data2[ariba.data2 == "yes"] <- 1
ariba.data3[ariba.data3 == "no"] <- 1
ariba.data3[ariba.data3 == "yes"] <- 0
ariba.data4[ariba.data4 == "no"] <- 1
ariba.data4[ariba.data4 == "yes"] <- 0
# 
# colnames(ariba.data) <- NULL
# colnames(ariba.data2) <- NULL

dbs.clusters <- read.table("DBS.clusters.tsv", header = F, sep = "\t")
pg.clusters <- read.table("PG.clusters.tsv", header = F, sep = "\t")

pg.clusters <- pg.clusters[pg.clusters$V1 %in% colnames(ariba.data4),]


pg.counts <- apply(ariba.data4,2, function(x) length(which(x != "NA")))

all.data <- cbind(ariba.data3, ariba.data4)

# gifsy.data <- cbind(Gifsys_DBS_snps[,c(2,3,4,5)], all.data[!rownames(all.data) %in% c("24","48","49"), ]) 


dd <- short.sample.data[,ncol(short.sample.data):1]

baps.data <- read.table("baps_metadata.txt",
                        sep = "\t", header = TRUE)

baps.3.lev <- read.table("../baps_data_3lev.txt",
                         sep = "\t", header = TRUE)
colnames(baps.3.lev) <- c("V2", "row", "layer", "cluster1", "cluster2", "cluster3")

baps.data <- merge(baps.data, baps.3.lev, by.x = "sample.y", by.y = "V2")

dd2 <- merge(dd, baps.data, by.x = "V2", by.y = "V2")

big.data <- merge(dd2, Gifsys_DBS_snps, by.x = "V2", by.y = "ID")

test.1 <- as.data.frame(lapply(ariba.data, as.numeric))
test.2 <- as.data.frame(lapply(ariba.data2, as.numeric))
test.3 <- as.data.frame(lapply(ariba.data3, as.numeric))
test.4 <- as.data.frame(lapply(ariba.data4, as.numeric))

test.2 <- test.2[,colSums(test.2) != 134]

resfinder.number <- as.data.frame(cbind(test.1[,1],rowSums(test.1[2:dim(test.1)[2]])))
vfdb.number <- as.data.frame(cbind(test.2[,1],rowSums(test.2[2:dim(test.2)[2]])))
large.dbs.number <- as.data.frame(cbind(test.3[,1],rowSums(test.3[2:dim(test.3)[2]])))
pg.number <- as.data.frame(cbind(test.4[,1],rowSums(test.4[2:dim(test.4)[2]])))

colnames(large.dbs.number) <- c("V2","no.DBS")
colnames(pg.number) <- c("V2","no.PG")
colnames(resfinder.number) <- c("V2","no.res")
colnames(vfdb.number) <- c("V2","no.vf")

big.data <- merge(big.data, large.dbs.number, by.x = "V2", by.y = "V2")
big.data <- merge(big.data, pg.number, by.x = "V2", by.y = "V2")
big.data <- merge(big.data, resfinder.number, by.x = "V2", by.y = "V2")
big.data <- merge(big.data, vfdb.number, by.x = "V2", by.y = "V2")
big.data <- merge(big.data, srst2.data, by.x = "V2", by.y = "Strain", all = T)


library(reshape2)
library(ggplot2)

new.data <- big.data[,c(23:35)]

new.data$cluster3 <- as.character(new.data$cluster3)
new.data$cluster1 <- as.character(new.data$cluster1)

# new.data[,5] <- round(as.numeric(new.data[,5]), 3)
# new.data[,4] <- round(as.numeric(new.data[,4]), 3)
# # new.data$cluster1[new.data$cluster1 == "epidemic"] <- "anthropogenic"
# 
# # new.data <- new.data[,-5]
# 
new.data <- new.data[,-c(1,2)]
new.data <- new.data[,-c(10,12)]
# # 
# # new.data$DBS <- round(as.numeric(as.character(new.data$DBS)),2)
# 
new.data$DBS <- as.numeric(as.character(new.data$DBS))
# # 
new.data <- new.data[!is.na(new.data$cluster3),]

melt_A <- melt(new.data)

library(ggpubr)
library(ggsignif)
melt_A$cluster3 <- factor(melt_A$cluster3, levels = c("1","2","3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19"))
# melt_A$cluster1 <- factor(melt_A$cluster1, levels = c("1","2"))


p1 <- ggplot(melt_A, aes(x=cluster3,y=value,fill=cluster3)) + geom_boxplot() + facet_wrap(~variable, scales = "free") #+
#   geom_signif(comparisons = list(c("1","2")),
#               map_signif_level=TRUE, test = "t.test")

# p2 <- ggplot(melt_A, aes(x=cluster1,y=value,fill=cluster1)) + geom_boxplot() + facet_wrap(~variable, scales = "free")
#
# ggsave("plot.2nd_group.pdf",plot = p1,  width = 600, height = 400, units = "mm",
#        dpi = 600, bg = "transparent")

# ggsave("plot.3rd_group_coreINV_core_genome.pdf",plot = p1,  width = 600, height = 400, units = "mm",
#        dpi = 600, bg = "transparent")

