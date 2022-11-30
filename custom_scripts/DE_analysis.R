#! /usr/bin/env R 
# written by Matthias Schmal
# distributed under GNU GPLv3

library("edgeR")

counts <- read.delim("count_combined.tsv")
design <- read.delim("design.tsv", header=FALSE, col.names=c("name", "group"))
hits <- read.delim("hits.txt", col.names="gene")

# sort design lexicographically
design <- design[order(design$name),]
g_counts <- as.data.frame(table(design$group))

# generate DGEList object containing the counts and group information
x <- DGEList(counts=counts[c(2:dim(counts)[2])], group=design$group, genes=counts[1])
# filter low-expressed transcripts, as they can not be used for meaningful statistics
keep <- filterByExpr(x)
x <- x[keep,,keep.lib.sizes=FALSE]

# DGE-analysis with replicates possible
if (1 %in% g_counts$Freq) {
    # calculate RPKM for each gene and transform to log2
    y <- cpm(x, log=TRUE)
    # calculate mean of each gene corresponding to their group tag
    y <- t(aggregate(t(y)[,2:dim(y)[1]], list(design$group), mean))


    # calculate foldchange for each gene
    foldchange <- y[,2] - y[,1]
    foldchange <- cbind(data.frame(x$genes[,]), data.frame(foldchange))
    foldchange <- foldchange[order(abs(foldchange$foldchange), decreasing=TRUE), ]
    colnames(foldchange) <- c("gene", "log2-foldchange")
    filtered <- foldchange[foldchange$gene %in% hits$gene, ]

    # write foldchanges as csv
    write.csv(foldchange, file="allFoldchange_sorted.csv", row.names=FALSE)
    write.csv(filtered, file="filteredFoldchange_sorted.csv", row.names=FALSE)
} else {
    # calculate scaling factors for calculation of effective library size
    x <- calcNormFactors(x)

    # read 2.10 in edgeR UserGuide
    x <- estimateCommonDisp(x)
    x <- estimateTagwiseDisp(x)

    # do test for DEG
    et <- exactTest(x)

    et <- cbind(x$genes, et$table)
    # sort according to PValue
    et <- et[order(abs(et$PValue)),]
    filtered <- et[et$genes %in% hits$gene, ]

    write.csv(et, file="allTags_sorted.csv", row.names=FALSE)
    write.csv(filtered, file="filteredTags_sorted.csv", row.names=FALSE)
}



