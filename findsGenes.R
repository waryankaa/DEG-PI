library(GEOquery)
library(limma)
library(umap)
library(tidyverse)

name <- read.delim("/Users/abbywaryanka/seniorComp/Streamlit/keyword.txt", header = FALSE)

# load series and platform data from GEO

name <- read.delim("/Users/abbywaryanka/seniorComp/Streamlit/keyword.txt", header = FALSE)

# load series and platform data from GEO
gset <- getGEO(name, GSEMatrix =TRUE, AnnotGPL=TRUE)
ending <- "_series_matrix.txt.gz"
type <- paste(name,ending)
type <- str_replace_all(string=type, pattern=" ", repl="")
platform <- gset[[type]]@annotation
if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# categorize samples
# sample_type_table <- select(gset@phenoData@data,title,geo_accession)

#read in files
sample_membership <- read.delim("/Users/abbywaryanka/seniorComp/Streamlit/group_membership.txt", header = FALSE)
sample_membership1 <- substring(sample_membership,2,)


# group membership for all samples
gsms <- paste0(sample_membership1)
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("normal", "cancer"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=50000)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05)

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))



write.csv(tT, "/Users/abbywaryanka/seniorComp/Streamlit/r.csv", row.names=FALSE)

