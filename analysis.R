## test

library(maftools)

## load inhibitor data
drugDat = read.csv("CCLE_NP24.2009_Drug_data_2015.02.24.csv", header = T)
vemDat = drugDat[drugDat[,"Compound"]=="PLX4720",]

## load mutation Data
mut = read.csv("CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf", sep = "\t", header =T, fill=T)
V600Mut = mut[grep("p.V600E", mut[,"Protein_Change"]), ]
BRAFMuts = V600Mut[V600Mut[,"Hugo_Symbol"]=="BRAF",]
skinBRAFMuts = as.character(BRAFMuts[grep("_SKIN", BRAFMuts[,"Tumor_Sample_Barcode"]), "Tumor_Sample_Barcode"])
  
## load expression data

read.gct <-
  function (file) {
    expr = read.table(file, skip = 2, header = TRUE, sep = "\t", quote = "")
    rownames(expr) = expr[,1]
    
    checkName = table(expr[,1])
    if(max(checkName) > 1) {
      stop(paste("Genes in gct file should be unique: ", names(which.max(checkName)), sep = " "))
    }
    expr = expr[,-c(1,2)]
    expr = as.matrix(expr)
    
    return(expr)
  }

expr = read.gct("CCLE_RNAseq_081117.rpkm.gct")

##
matchedDrug = vemDat[match(skinBRAFMuts, vemDat[,"CCLE.Cell.Line.Name"]),]
matchedDrug = matchedDrug[!is.na(matchedDrug[,1]),]

resistThresh = 5
sensThresh = 0.5

resistantCL = matchedDrug[matchedDrug[,"IC50..uM."]>=resistThresh, "CCLE.Cell.Line.Name"]
sensitiveCL = matchedDrug[matchedDrug[,"IC50..uM."]<=sensThresh, "CCLE.Cell.Line.Name"]

## narrow down expression matrix

## strip ENSMBLE ID
rownames(expr) = sub("\\..*", "", rownames(expr))

resExpr = expr[, colnames(expr) %in% resistantCL]
senExpr = expr[, colnames(expr) %in% sensitiveCL]

myGenes = read.csv("top_genes.csv", header=T, sep="\t")
myGenes = unique(myGenes[, "gene"])

library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset = "hsapiens_gene_ensembl")
ensmblGenes = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "hgnc_symbol", values = myGenes, mart = ensembl)
ensmblGenes = as.data.frame(ensmblGenes)

## generate_boxplots

for (i in 1:nrow(ensmblGenes)){
    
    currGeneEn =  ensmblGenes[i,1] 
    currGeneHg = ensmblGenes[i,2]
    
    ## get expression
    currResExp = resExpr[currGeneEn,]
    currSenExp = senExpr[currGeneEn,]
    
    currExp = cbind(as.data.frame(c(currResExp, currSenExp)), Type=factor(c(rep("Resistant", length(currResExp)), rep("Sensitive", length(currSenExp) ))) ) 
    colnames(currExp)[1] = "Exp" 
      
    tiff(filename=paste("./figures/", currGeneHg, ".tiff", sep=""), width = 5, height = 5, units = 'in', res = 250)
    ## make boxplot
    boxplot(Exp~Type,data=currExp, lwd=2, col=c('powderblue', 'mistyrose'), outline=F, ylab="RPKM", main=currGeneHg)
    currExpnA375 = currExp[-which(rownames(currExp)=="A375_SKIN"),]
    stripchart(Exp ~ Type, vertical = TRUE, data = currExpnA375, 
               method = "jitter", add = TRUE, pch = 20, col = 'royalblue', cex=2)
    text(x=2, y=currExp["A375_SKIN",1], labels = "A375")
    dev.off()
    
}

## generate_scatter
allDat = c()
for (i in 1:nrow(ensmblGenes)){
  
  currGeneEn =  ensmblGenes[i,1] 
  currGeneHg = ensmblGenes[i,2]
  
  ## get expression
  IC50s = matchedDrug[,c(1,11)]
  mcExpr = expr[currGeneEn, colnames(expr) %in% IC50s[,1]]
  
  finalExpr = mcExpr[match(IC50s[,1], names(mcExpr))]
  
  allDat = rbind(allDat, cbind( as.data.frame(cbind(IC50s, finalExpr)), type=factor(rep(currGeneHg, nrow(IC50s))) ) )
  
}

colnames(allDat) = c("name", "IC50", "RPKM","type")

p = ggplot(allDat, aes(IC50, RPKM))
p+geom_point(aes(colour=factor(allDat[,"type"])), size=2)+ylim(0,100)


