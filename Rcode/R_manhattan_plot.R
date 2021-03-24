library(ggplot2)
library(cowplot)

# read in fst file in local computer
fstFull = read.delim("SN.fst")
fstSub <- fstFull[!is.na(fstFull$WEIR_AND_COCKERHAM_FST),]

# prepare for ploting each chromosome
chr01p <- ggplot(fstSub[fstSub$CHROM=="Chr01",], aes(x=POS, y=WEIR_AND_COCKERHAM_FST)) + geom_point() + labs(x="Chr 01 Pos", y="Fst") + ylim(0,0.7)
chr02p <- ggplot(fstSub[fstSub$CHROM=="Chr02",], aes(x=POS, y=WEIR_AND_COCKERHAM_FST)) + geom_point() + labs(x="Chr 02 Pos", y="Fst") + ylim(0,0.7)
chr03p <- ggplot(fstSub[fstSub$CHROM=="Chr03",], aes(x=POS, y=WEIR_AND_COCKERHAM_FST)) + geom_point() + labs(x="Chr 03 Pos", y="Fst") + ylim(0,0.7)
chr04p <- ggplot(fstSub[fstSub$CHROM=="Chr04",], aes(x=POS, y=WEIR_AND_COCKERHAM_FST)) + geom_point() + labs(x="Chr 04 Pos", y="Fst") + ylim(0,0.7)
chr05p <- ggplot(fstSub[fstSub$CHROM=="Chr05",], aes(x=POS, y=WEIR_AND_COCKERHAM_FST)) + geom_point() + labs(x="Chr 05 Pos", y="Fst") + ylim(0,0.7)
chr06p <- ggplot(fstSub[fstSub$CHROM=="Chr06",], aes(x=POS, y=WEIR_AND_COCKERHAM_FST)) + geom_point() + labs(x="Chr 06 Pos", y="Fst") + ylim(0,0.7)
chr07p <- ggplot(fstSub[fstSub$CHROM=="Chr07",], aes(x=POS, y=WEIR_AND_COCKERHAM_FST)) + geom_point() + labs(x="Chr 07 Pos", y="Fst") + ylim(0,0.7)
chr08p <- ggplot(fstSub[fstSub$CHROM=="Chr08",], aes(x=POS, y=WEIR_AND_COCKERHAM_FST)) + geom_point() + labs(x="Chr 08 Pos", y="Fst") + ylim(0,0.7)
chr09p <- ggplot(fstSub[fstSub$CHROM=="Chr09",], aes(x=POS, y=WEIR_AND_COCKERHAM_FST)) + geom_point() + labs(x="Chr 09 Pos", y="Fst") + ylim(0,0.7)
chr10p <- ggplot(fstSub[fstSub$CHROM=="Chr10",], aes(x=POS, y=WEIR_AND_COCKERHAM_FST)) + geom_point() + labs(x="Chr 10 Pos", y="Fst") + ylim(0,0.7)
chr11p <- ggplot(fstSub[fstSub$CHROM=="Chr11",], aes(x=POS, y=WEIR_AND_COCKERHAM_FST)) + geom_point() + labs(x="Chr 11 Pos", y="Fst") + ylim(0,0.7)
chr12p <- ggplot(fstSub[fstSub$CHROM=="Chr12",], aes(x=POS, y=WEIR_AND_COCKERHAM_FST)) + geom_point() + labs(x="Chr 12 Pos", y="Fst") + ylim(0,0.7)
chr13p <- ggplot(fstSub[fstSub$CHROM=="Chr13",], aes(x=POS, y=WEIR_AND_COCKERHAM_FST)) + geom_point() + labs(x="Chr 13 Pos", y="Fst") + ylim(0,0.7)
chr14p <- ggplot(fstSub[fstSub$CHROM=="Chr14",], aes(x=POS, y=WEIR_AND_COCKERHAM_FST)) + geom_point() + labs(x="Chr 14 Pos", y="Fst") + ylim(0,0.7)
chr15p <- ggplot(fstSub[fstSub$CHROM=="Chr15W",], aes(x=POS, y=WEIR_AND_COCKERHAM_FST)) + geom_point() + labs(x="Chr 15W Pos", y="Fst") + ylim(0,0.7)
chr16p <- ggplot(fstSub[fstSub$CHROM=="Chr16",], aes(x=POS, y=WEIR_AND_COCKERHAM_FST)) + geom_point() + labs(x="Chr 16 Pos", y="Fst") + ylim(0,0.7)
chr17p <- ggplot(fstSub[fstSub$CHROM=="Chr17",], aes(x=POS, y=WEIR_AND_COCKERHAM_FST)) + geom_point() + labs(x="Chr 17 Pos", y="Fst") + ylim(0,0.7)
chr18p <- ggplot(fstSub[fstSub$CHROM=="Chr18",], aes(x=POS, y=WEIR_AND_COCKERHAM_FST)) + geom_point() + labs(x="Chr 18 Pos", y="Fst") + ylim(0,0.7)
chr19p <- ggplot(fstSub[fstSub$CHROM=="Chr19",], aes(x=POS, y=WEIR_AND_COCKERHAM_FST)) + geom_point() + labs(x="Chr 19 Pos", y="Fst") + ylim(0,0.7)

# output figure to local computer
fullFigure <- plot_grid(chr01p,chr02p,chr03p,chr04p,chr05p,chr06p,chr07p,chr08p,chr09p,chr10p,chr11p,chr12p,chr13p,chr14p,chr15p,chr16p,chr17p,chr18p,chr19p,align="hv",ncol=4)
ggsave(filename="fst.png", height=30, width=18)

# high-Fst locus selection
  # Note: below the number '0.995' can be manipulated if required. It is the quantile of all Fst values among genome.
  # If your output is too many and not clean (spreaded in several chromosomes/over 150 loci), you can choose higher quantile
  # If your output is lower than 30 loci you may lower the quantile to 0.99 or even lower.
outFrame <- fstSub[fstSub$WEIR_AND_COCKERHAM_FST > quantile(fstSub$WEIR_AND_COCKERHAM_FST, 0.995), c(1,2)]
write.table(outFrame, 'includeCoords.txt', sep="\t", eol="\n", quote=F, row.names=F, col.names=F)
