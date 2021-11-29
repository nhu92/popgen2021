# MK test
# adapted & edited from: https://eacooper400.github.io/gen8900/exercises/mk.html
# Nan Hu
# Mar. 27, 2021

# Install packages from BiocManager/BiocLite
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")
BiocManager::install("VariantAnnotation")
BiocManager::install("GenomicFeatures")
BiocManager::install("Rsamtools")
library(dplyr)
library(devtools)
# load some functions from the author
source_gist("c954c757db02eb6c4ccfae1aa090658c", filename="compGen17_fxns.R", quiet=TRUE)

# read input file
library(Rsamtools)
fa = open(FaFile("Salix_purpurea_var_94006.mainGenome.fasta"))

library(GenomicFeatures)
txdb=makeTxDbFromGFF(file="Spurpurea94006v5.1.gene_exons_NoZ.gff3", format="gff3")

library(VariantAnnotation)
vcf.object=readVcf(file="sreticulata.vcf")

# predict coding positions
effects = predictCoding(vcf.object, txdb, fa)

# match up gene names
id.from.effects = names(ranges(effects))
id.from.vcf = rownames(info(vcf.object))
m = match(id.from.vcf, id.from.effects)

info(vcf.object)$CSQ = effects$CONSEQUENCE[m]
# names(vcf.object) = effects$GENEID[m]

ij <- effects$GENEID[m] 
ij[is.na(ij)] <- 'nongeneric'
names(vcf.object) <- ij


writeVcf(vcf.object, "sreticulata_annot.vcf")

# re-write a function that read vcf files
read.vcf_m <- function(file, special.char="##", ...) {
  my.search.term=paste0(special.char, ".*")
  all.lines=readLines(file)
  clean.lines=gsub(my.search.term, "",  all.lines)
  clean.lines=gsub("#CHROM", "CHROM", clean.lines)
  read.table(..., text=paste(clean.lines, collapse="\n"), fill = TRUE)
}

# some data cleaning and filtering
my.data=read.vcf_m("sreticulata_annot.vcf", header=TRUE, stringsAsFactors=FALSE)
head(my.data)
my.data = my.data[which(my.data$ID != 'nongeneric'),]
my.data = my.data[!(is.na(my.data$POS)),]

head(my.data)
my.data$INFO = gsub(".*;CSQ=", "", my.data$INFO)
my.data = my.data[which(my.data$INFO != 'not'),]
head(my.data)

# new table for calculating reference allele frequencies
new <- my.data[,c(3,8)]
head(new)
new$refAF <- NA
for (i in 1:nrow(my.data)) {
  my.row <- as.vector(my.data[i,], mode="character")
  genotypes <- get.field(my.row[10:length(my.row)], my.row[9], "GT")
  all.counts <- count.genotypes(genotypes)
  ref.freq <- allele.freq(all.counts)["p"]
  new$refAF[i] <- ref.freq
}
# based on allele frequencies of reference allele, determine site classes
new$SiteClass <- new$refAF
new$SiteClass[which(new$refAF==0)] <- "Fixed"
new$SiteClass[which(new$refAF>0)] <- "Polymorphic"
new <- new[which(new$refAF != 1), ]
siteclass.tab <- new[,c(1,2,4)]

# unique gene list
geneList=unique(my.data$ID)

# result table format
my.results = data.frame(Genes=geneList, Fn=rep(0, length(geneList)), Fs=rep(0, length(geneList)), Pn=rep(0, length(geneList)), Ps=rep(0, length(geneList)))

# calc Fn Fs Pn Ps
for (i in 1:length(geneList)) {
  gene=siteclass.tab[which(siteclass.tab$ID==geneList[i]),]
  syn=which(gene$INFO=="synonymous")
  non=which(gene$INFO!="synonymous")
  my.results$Fn[i] = length(intersect(non, which(gene$SiteClass=="Fixed")))
  my.results$Fs[i] = length(which(syn %in% which(gene$SiteClass=="Fixed")))
  my.results$Pn[i] = nrow(subset(gene, (gene$SiteClass=="Polymorphic" & gene$INFO!="synonymous")))
  my.results$Ps[i] = nrow(gene) - sum(my.results[i,2:5])
}

# MK test
my.results$BetweenN_S = my.results$Fn/my.results$Fs
my.results$WithinN_S = my.results$Pn/my.results$Ps
my.results$NI = (my.results$Pn/my.results$Ps)/(my.results$Fn/my.results$Fs)
my.results$pValues=apply(my.results[,2:5], 1, function(x) fisher.test(matrix(x, ncol=2))$p.value)

head(my.results,15)
# save into file
write.table(my.results, 
            file = "sreticulata.mktest.txt",
            sep = "\t", 
            row.names = FALSE,
            col.names = TRUE)
