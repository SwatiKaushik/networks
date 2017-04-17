# Program to calculate the network biomarkers of a given kinase in patient's omics data and
# validate if the networks can predict drug response in PDX models
# Swati Kaushik
######################################### OUTPUT FILES ###################################

# unique interactors OUT1: "unique.kinase.interactors"
# gene expression of interactors OUT2: "Gene.expression.interactors"
# samples with mutations in kinase OUT3: "mutated.sample.set"
# t.test/fdr results OUT4: "t.test.results"
# network activity data OUT5: "Network-activity-data"
# Heatmap of network activity OUT6: "Network-activity-heatmap.pdf"

##########################################################################################

rm(list=ls())
library(gplots)
library(devtools)
library(plotly)

setwd("/users/swati/Desktop/NetworkActivity")

#function calculates network activity of the selected gene

network.activity.calculations <- function (x, p.count){
	output <- vector()
	p <- as.vector(colnames(x))

	for (i in 1:ncol(x)){
		sum.cols <- sum(x[1:p.count,i])
		sub.cols <- sum(x[(p.count+1):nrow(x),i])
		network.activity <- sum.cols - sub.cols
		no.data <- paste(p[i], network.activity,collapse='	')
		output <- rbind(output, data.frame(p[i],network.activity))
	}
	rownames(output) <- output$p.i.
	write.table(output[,2], file = "NA.out", sep="\t", row.names = rownames(output))
	return (output)
}

#calculates z score of the given vector
z.score <- function(x){
	na.vector <- t(x)
	na.vector.t <- scale(t(na.vector))
	attr(na.vector.t, "scaled:scale") = NULL
	attr(na.vector.t, "scaled:center") = NULL
	na.vector <- round(t(na.vector.t), digits=3)
	return (na.vector)
}

NA.component.matrix <- function(na.vector1, pn.matrix){
	colnames(na.vector1) <- colnames(pn.matrix)
	new.matrix <- rbind(pn.matrix, na.vector1)
	names.matrix <- rownames(pn.matrix)
	appended.names.matrix <- append(names.matrix, "EGFR-NA")
	rownames(new.matrix) <- appended.names.matrix
	transposed.new.matrix <- as.data.frame(t(new.matrix))
	transposed.new.matrix <- transposed.new.matrix[order(transposed.new.matrix$`EGFR-NA`),]
	transposed.new.matrix.sorted <- as.data.frame(t(transposed.new.matrix))
	return(transposed.new.matrix.sorted)
}

###################################### INPUT FILES #######################################

#name of the gene for which network activity has to be calculated
kinase.name <- 'EGFR'

### input files required
cat ("\nreading input file......\n")
#proteomics data
proteomics.file.inhouse <- read.table("proteomics-data-inhouse.txt", header=T, sep="\t")
proteomics.file.public <- read.table("proteomics-data-public.txt", header=T, sep="\t")

#cancer patients gene expression and mutation data
gene.expression <- read.table("final.set.IDname1.T.log.Zscore", header=T, sep="\t")
mutation.patients <- read.table("luad.mutation.maf", header=T, sep="\t")

#test on pdx data 
pdx.expression <- read.table("NSCLc.pdx.exp", header=T, sep="\t")
pdx.response.data <- read.table("pdx-response-data.txt", header=T, sep="\t")


######################################## Cancer patients datasets ########################

### process proteomics datasets
cat ("processing proteomics datasets\n")
###

# get proteomics data file (developed in house)
proteomics.kinase.data <- proteomics.file.inhouse[grep(kinase.name,proteomics.file.inhouse$Bait ),]

#put a cutoff to get the interactors of the kinase
proteomics.subset.zs <- subset(proteomics.kinase.data$Prey, proteomics.kinase.data$ZscoreWD >2.5)
unq.proteomics.subset.zs <- levels(droplevels(proteomics.subset.zs))

# open proteomics data file (publicly available)
proteomics.kinase.public.data <- proteomics.file.public[grep(kinase.name,proteomics.file.public$Bait ),]
unq.proteomics.set2 <- levels(droplevels(proteomics.kinase.public.data$Prey))

##get unique interactors of the kinase from the two files
unique.ints <- unique(c(unq.proteomics.set2, unq.proteomics.subset.zs))
unique.ints.kinase <- append(kinase.name, unique.ints)
unique.interactors <- paste(paste("^", unique.ints.kinase, sep=""), "$", sep="")
write.table(unique.ints.kinase, file="unique.kinase.interactors", sep="\t")

### process gene expression data of lung cancer patients
cat ("processing gene expression dataset\n")
###

rownames(gene.expression) <- sub('\\..*', '', rownames(gene.expression))
interactor.expression <- gene.expression[grep(paste(unique.interactors , collapse="|"), rownames(gene.expression)),]
colnames.substring <- substr(colnames(interactor.expression),1,15)
colnames.substring.2 <- gsub(".",'-',colnames.substring, fixed=T)
colnames(interactor.expression) <- colnames.substring.2
write.table(interactor.expression, file="Gene.expression.interactors", sep="\t")

### process mutation data of cancer patients
cat ("processing mutation dataset\n")
###

kinase.mutation <- mutation.patients[grep(kinase.name, mutation.patients$Hugo_Symbol),]
kinase.mutation.type <- kinase.mutation[grep("Missense_Mutation|In_Frame_Del", kinase.mutation$Variant_Classification),]
mutated.samples <- levels(droplevels(kinase.mutation.type$Tumor_Sample_Barcode))
mutated.samples.names <- substr(mutated.samples,1,15)

mutated.set <- interactor.expression[intersect(colnames(interactor.expression), mutated.samples.names)]
mutated.set <- mutated.set[order(rownames(mutated.set)),]
wild.set <- interactor.expression[setdiff(colnames(interactor.expression), mutated.samples.names)]
wild.set <- wild.set[order(rownames(wild.set)),]
write.table(mutated.set, file="mutated.sample.set", sep="\t")
#write.table(wild.set, file="wild.set", sep="\t")

### Identify interactors correlated with mutations of the kinase
cat ("calculating correlated interactors\n")
###for each interactor compare the expression of EGFR mutant set with wild type set

t.testresults <- vector("list", nrow(mutated.set))

for (j in seq(nrow(mutated.set))){
  t.testresults[[j]] <- t.test(mutated.set[j,], wild.set[j,])
}

get.t.test <- as.data.frame(t(sapply(t.testresults, function(x) unlist(x[c("estimate","p.value","statistic")]))))
rownames(get.t.test) <- rownames(mutated.set)
colnames(get.t.test) <- c("exp.mean-EGFR-mutant", "exp.mean-EGFR-wild","pvalue","statistics")
get.t.test$fdr <- p.adjust(get.t.test$pvalue, method="fdr")

#perform FDR
sorted.get.t.test <- get.t.test[order(get.t.test$fdr),]
write.table(sorted.get.t.test, file="t.test.results", sep="\t", quote=FALSE)
positive.corr <- subset(rownames(sorted.get.t.test) , sorted.get.t.test$`exp.mean-EGFR-mutant` >0 & sorted.get.t.test$fdr < 0.05)
negative.corr <- subset(rownames(sorted.get.t.test) , sorted.get.t.test$`exp.mean-EGFR-mutant` <0 & sorted.get.t.test$fdr < 0.05)
positive.expression <- interactor.expression[grep(paste(positive.corr, collapse = "|"), rownames(interactor.expression)),]
negative.expression <- interactor.expression[grep(paste(negative.corr, collapse = "|"), rownames(interactor.expression)),]

pn.matrix <- rbind(positive.expression, negative.expression)
order.corr <- c(positive.corr, negative.corr)
pn.matrix <- pn.matrix[match(order.corr, rownames(pn.matrix)),]

###calculate network activity using correlated interactors
cat ("calculating network activity\n")
###

p.count <- nrow(positive.expression)
na.output <- network.activity.calculations(pn.matrix, p.count)

# calculate z score of NA
na.vector1 <- z.score(na.output[,2])

#create all NA component matrix
transposed.new.matrix.sorted <- NA.component.matrix (na.vector1, pn.matrix)
write.table(transposed.new.matrix.sorted, file = "Network-activity-data", sep="\t")

#sort mutation vector according to transposed matrix
m.vector <- colnames(transposed.new.matrix.sorted) %in% colnames(mutated.set)
column_col1 <- ifelse(m.vector == "TRUE", "red", "gray65")
column_annotation1 <- as.matrix(column_col1)

###plot heatmap

#Load latest version of heatmap.3 function
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

#pdf(file="Network-activity-heatmap.pdf")

heatmap.3(as.matrix(transposed.new.matrix.sorted),
			col = bluered(256),
			ColSideColors = column_annotation1, 
			Rowv = FALSE, 
			Colv = FALSE, 
			dendrogram = "none",
			labCol = FALSE,
			margins = c(6,12),
			ColSideColorsSize = 1)
#dev.off()	


############### Check performance of network signature on PDX data ###########################

# get PDX data 
order.interactors <- paste(paste("^", order.corr, sep=""), "$", sep="")
pdx.subset.matrix <- pdx.expression[grep(paste(order.interactors , collapse="|"), rownames(pdx.expression)), ]
pdx.int.matrix <- pdx.subset.matrix[match(order.corr, rownames(pdx.subset.matrix)),]

#calculate network activity of PDX
pdx.na.output <- network.activity.calculations(pdx.int.matrix, 14)

# calculate z score of NA
na.vector1 <- z.score(pdx.na.output[,2])

#create all NA component matrix
transposed.new.matrix.sorted <- NA.component.matrix (na.vector1, pdx.int.matrix)
#write.table(transposed.new.matrix.sorted, file = "pdx.new.entire.matrix.out", sep="\t")

pdx.response.na <- merge(pdx.na.output, pdx.response.data, all=T, by='row.names')

#calculate correlation of network activity with drug response data

pdx.corr <- cor(pdx.response.na$network.activity, pdx.response.na$Bestavgresponse, method="spearman")
krasmutant <- c("X.1156", "X.1442", "X.2082","X.2094","X.4819", "X.1172", "X.1835")
kras.data <- pdx.response.na[grep(paste(krasmutant, collapse="|"), pdx.response.na$p.i.),]
egfrmutant <- c("X.1683")
EGFR.data <- pdx.response.na[grep(paste(egfrmutant, collapse="|"), pdx.response.na$p.i.),]

#add mutations to pdx data
kras.add <- rownames(pdx.response.na) %in% row.names(kras.data)
kras.add.list <- ifelse(kras.add == "TRUE", "kras", "none")
egfr.no <- which(rownames(pdx.response.na) %in% row.names(EGFR.data) == TRUE)
kras.add.list[egfr.no] <- "EGFR"
pdx.plot.data <- cbind(pdx.response.na, kras.add.list)
row <- pdx.plot.data$Row.names

# Plot PDX data
# scatter plot data 
pal <- c("red", "blue", "green")
y <- list(
  title = "EGFR Network activity"
)
x <- list(
  title = "Change in tumor volume when treated with drug (%)"
)
p <- plot_ly(
	data = pdx.plot.data, x = ~Bestavgresponse, y = ~network.activity, 
	text = ~paste("sample:", row),
	color = ~kras.add.list, colors = pal) %>%	
layout(xaxis =x, yaxis = y)

suppressWarnings(print(p))

#boxplot responders vs. nonresponders
sorted.pdx.plot.data <- pdx.plot.data[order(pdx.plot.data$Bestavgresponse),]
sorted.pdx.plot.data$response <- c(rep("responder",10), rep("nonresponder", 15))
sorted.pdx.plot.data$color <- c(rep("green",10), rep("yellow", 15))
p <- plot_ly(sorted.pdx.plot.data, x= ~response,y = ~network.activity, color = ~color, type = "box")
suppressWarnings(print(p))

pdx.pvalue <- wilcox.test(sorted.pdx.plot.data$network.activity[1:10], sorted.pdx.plot.data$network.activity[11:25])$p.value
pdx.pvalue 

