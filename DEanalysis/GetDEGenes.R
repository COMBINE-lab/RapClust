args <- commandArgs(trailingOnly=TRUE)
methods <- lapply(args[1], tolower)
quantdir <- args[2]
clustFile <- args[3]
outdir <- args[4]
species <- args[5]

library(iCOBRA)
suppressMessages(library(DESeq2))
suppressMessages(library(tximport))
suppressMessages(library(data.table))
suppressMessages(library(readr))

getTXImport <- function(methodName , species) {
    # Where the quantification results are
    quantDir <- normalizePath(file.path(quantdir))
    message("quantDir is ", quantDir) 

    # Process depending on the method
    tx2gene <- read.csv(file.path(clustFile), sep="\t")
    if (methodName == "sailfish")
	tx2gene <- tx2gene[,c(2,1)]
    message(clustFile)
    message("import sailfish results")
    files <- file.path(quantDir, dir(file.path(quantDir)), "quant.sf")
    message(files)
    if (species == "chicken") {
        names(files) <- paste0("sample", c("A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4")) 
    } else {
        names(files) <- paste0("sample", c("A1", "A2", "A3", "B1", "B2", "B3")) 
    }
    message(names(files))
    txi <- tximport(files, type="sailfish", countsFromAbundance = "scaledTPM", tx2gene = tx2gene)
    
    return(txi)
} 

getDEGenesFromTxi <- function(txi, meta, cond_name, level1, level2, sample_name, species) {
    if (species == "chicken") {    
        sampleTable <- data.frame(condition=factor(rep(c("A","B"),each=4)))
    } else {
           sampleTable <- data.frame(condition=factor(rep(c("A","B"),each=3)))
    }
    rownames(sampleTable) <- colnames(txi$counts)
    message(sampleTable)
	
    txi$counts <- round(txi$counts)
    keep_feat <- rownames(txi$counts[rowSums(is.na(txi$counts)) == 0 & rowSums(txi$counts) != 0, ])
    txi <- lapply(txi, function(w) {
      if (!is.null(dim(w))) w[match(keep_feat, rownames(w)), ]
      else w
      })
    colData = meta[match(colnames(txi$counts), meta[, sample_name]), ]
    dsd <- DESeq2::DESeqDataSetFromTximport(txi, 
                                    colData = colData, 
                                    design = as.formula(paste0("~", cond_name)))
   
    dsd <- DESeq(dsd, test="Wald", fitType="local", betaPrior=TRUE)
    res <- as.data.frame(results(dsd, contrast = c("condition", "B", "A"), cooksCutoff = FALSE, independentFiltering = FALSE))
    #setorder(res, -padj)
    return(res)
}

if (species == "chicken") {
    meta <- data.frame(sample = paste0("sample", c("A1", "A2", "A3","A4", "B1", "B2", "B3", "B4")),
                       condition = c("A", "A", "A", "A", "B", "B", "B", "B"),
                       stringsAsFactors = FALSE)
} else {
    meta <- data.frame(sample = paste0("sample", c("A1", "A2", "A3", "B1", "B2", "B3")),
                       condition = c("A", "A", "A", "B", "B", "B"),
                       stringsAsFactors = FALSE)
}
rownames(meta) <- meta$sample

first = TRUE
for (m in methods) {

  if ((m == "sailfish") | (m == "truth") || (m == "cdhit")) {
	txis <- lapply(methods, getTXImport, species)
	names(txis) <- methods

	allres <- lapply(txis, getDEGenesFromTxi, meta = meta, cond_name = "condition", level1 = "A", level2 = "B", sample_name = "sample", species)
	names(allres) <- methods

	options("scipen"=100, "digits"=4)
    	padj <- data.frame(sailfish = allres[[m]]$padj, row.names = rownames(allres[[m]]))
	format(padj, scientific=FALSE)
	write.table(padj, file = file.path(outdir, paste0(m,"padj.txt")), sep = "\t", row.names = rownames(allres[[m]]), quote = FALSE, col.names = FALSE)
  } else if (m == "corset") {
	dat <- data.matrix(read.csv(file.path(clustFile), header = TRUE, row.names = 1, sep = "\t"))
	clusters <- rownames(dat)	
	rownames(dat) <- NULL
	
    if (species == "chicken") {
	    condition=factor(rep(c("A","B"),each=4))
    } else {
        condition=factor(rep(c("A","B"),each=3))
    }
    	dsd <- DESeq2::DESeqDataSetFromMatrix(dat, DataFrame(condition), ~ condition)
  
    	dsd <- DESeq(dsd, test="Wald", fitType="local", betaPrior=TRUE)
	
    	allres <- as.data.frame(results(dsd, contrast = c("condition", "B", "A"), cooksCutoff = FALSE, independentFiltering = FALSE))
	#names(allres) <- methods

    	padj <- data.frame(corset = allres$padj, row.names = clusters)
	format(padj, scientific=FALSE)
	write.table(padj, file = file.path(outdir, paste0(m, "padj.txt")), sep = "\t", row.names = clusters, quote = FALSE, col.names = FALSE)
    	message("Loading corset data")
  }
  first = FALSE
}

