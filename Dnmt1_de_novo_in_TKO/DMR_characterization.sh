# Methylation distribution DMRs and CRs
for i in WT.bw TKOlike.bw DKOzero_P1.bw DKOzero_P5.bw DKOzero_P15.bw DKOzero_P25.bw TKO.bw TKO_Dnmt1PiggyBac.bw TKO_Dnmt1_pcat.bw; do
    out=`basename ${i} | sed 's/.bw//'`
    echo ${i}
    bigWigAverageOverBed ${i} ../DKO_TKO_DMRs_numerated_v3.bed ${out}_DMR.txt
    bigWigAverageOverBed ${i} ../DKO_TKO_CRs_numerated_v3.bed ${out}_CR.txt
done

```R
require(plyr)
require(ggplot2)
require(ggridges)
require(reshape2)
require(vioplot)

mESC_DMR <- read.table('WT_DMR.txt', header=F)
TKOlike_DMR <- read.table('TKOlike_DMR.txt', header=F)
DKOzero_P1_DMR <- read.table('DKOzero_P1_DMR.txt', header=F)
DKOzero_P5_DMR <- read.table('DKOzero_P5_DMR.txt', header=F)
DKOzero_P15_DMR <- read.table('DKOzero_P15_DMR.txt', header=F)
DKOzero_P25_DMR <- read.table('DKOzero_P25_DMR.txt', header=F)
TKO_DMR <- read.table('TKO_DMR.txt', header=F)
TKO_Dnmt1PiggyBac_DMR <- read.table('TKO_Dnmt1PiggyBac_DMR.txt', header=F)
TKO_Dnmt1_pcat_DMR <- read.table('TKO_Dnmt1_pcat_DMR.txt', header=F)
DMR <- Reduce(function(x, y) merge(x, y, all.x=T, all.y=T, by=c('V1')), list(mESC_DMR[,c(1,6)], TKOlike_DMR[,c(1,6)], DKOzero_P1_DMR[,c(1,6)], DKOzero_P5_DMR[,c(1,6)], DKOzero_P15_DMR[,c(1,6)], DKOzero_P25_DMR[,c(1,6)], TKO_DMR[,c(1,6)], TKO_Dnmt1PiggyBac_DMR[,c(1,6)], TKO_Dnmt1_pcat_DMR[,c(1,6)]))

mESC_CR <- read.table('WT_CR.txt', header=F)
TKOlike_CR <- read.table('TKOlike_CR.txt', header=F)
DKOzero_P1_CR <- read.table('DKOzero_P1_CR.txt', header=F)
DKOzero_P5_CR <- read.table('DKOzero_P5_CR.txt', header=F)
DKOzero_P15_CR <- read.table('DKOzero_P15_CR.txt', header=F)
DKOzero_P25_CR <- read.table('DKOzero_P25_CR.txt', header=F)
TKO_CR <- read.table('TKO_CR.txt', header=F)
TKO_Dnmt1PiggyBac_CR <- read.table('TKO_Dnmt1PiggyBac_CR.txt', header=F)
TKO_Dnmt1_pcat_CR <- read.table('TKO_Dnmt1_pcat_CR.txt', header=F)
CR <- Reduce(function(x, y) merge(x, y, all.x=T, all.y=T, by=c('V1')), list(mESC_CR[,c(1,6)], TKOlike_CR[,c(1,6)], DKOzero_P1_CR[,c(1,6)], DKOzero_P5_CR[,c(1,6)], DKOzero_P15_CR[,c(1,6)], DKOzero_P25_CR[,c(1,6)], TKO_CR[,c(1,6)], TKO_Dnmt1PiggyBac_CR[,c(1,6)], TKO_Dnmt1_pcat_CR[,c(1,6)]))colnames(CR) <- c('ID','CR_mESC','CR_TKOlike','CR_DKOzero_P1','CR_DKOzero_P5','CR_DKOzero_P15','CR_DKOzero_P25','CR_KH2_TKO_C5_P35','CR_KH2_TKO_C5_p43','CR_KH2_TKO_Dnmt1PiggyBac_C37_P34_d17','CR_KH2_TKO_Dnmt1PiggyBac_C37_p57','CR_KH2_TKO_Dnmt1_pcat_C44_p40')

colnames(DMR) <- gsub('_DMR','',colnames(DMR))
colnames(CR) <- gsub('_CR','',colnames(CR))

pdf('figures/DMR_CR_splitviolin.pdf', width=21)
vioplot(DMR[,-1], cex.axis=0.25, side = "left", col='yellow')
vioplot(CR[,-1], cex.axis=0.25, side = "right", add = T, col='grey')
dev.off()
```


# Average methylation profile plot + heatmap
## DMRs
deeptools/bin/computeMatrix scale-regions \
-S WT.bw DKOzero_P15.bw TKOlike.bw \
-R DKO_TKO_DMRs_numerated_v1.bed \
-o DKO_TKO_DMRs.tab.gz \
-a 5000 -b 5000 \
--binSize 50

deeptools/bin/plotProfile \
-m DKO_TKO_DMRs.tab.gz \
-out figures/DKO_TKO_DMRs_numerated_v1_profile.pdf \
--plotType std \
--yMin -0.1 \
--yMax 1.1 \
--startLabel DMR \
--endLabel DMR \
--perGroup \
--colors 'grey' 'orange' 'green'

## CRs
deeptools/bin/computeMatrix scale-regions \
-S WT.bw DKOzero_P15.bw TKOlike.b \
-R DKO_TKO_CRs_numerated_v1.bed \
-o DKO_TKO_CRs.tab.gz \
-a 5000 -b 5000 \
--binSize 50

deeptools/bin/plotProfile \
-m DKO_TKO_CRs.tab.gz \
-out figures/DKO_TKO_CRs_numerated_v1_profile.pdf \
--plotType std \
--yMin -0.1 \
--yMax 1.1 \
--startLabel CR \
--endLabel CR \
--perGroup \
--colors 'grey' 'orange' 'green'


# Clonal behaviour
bedtools unionbedg -header -names scIAP-1 scIAP-2 scIAP-3 scIAP-4 scIAP-5 scIAP-6 scIAP-7 scIAP-8 scIAP-9 -i DKOzero-p5-scIAP-1.bedgraph DKOzero-p5-scIAP-2.bedgraph DKOzero-p5-scIAP-3.bedgraph DKOzero-p5-scIAP-4.bedgraph DKOzero-p5-scIAP-5.bedgraph DKOzero-p5-scIAP-6.bedgraph DKOzero-p5-scIAP-7.bedgraph DKOzero-p5-scIAP-8.bedgraph DKOzero-p5-scIAP-9.bedgraph >clones.bedgraph

```R
require(ggplot2)
require(reshape2)

data <- read.table('clones.bedgraph', header=T)

pdf('figures/clones.pdf')
ggplot(melt(data[,-2], id.vars=c('chrom','end')), aes(x=end, y=value, color=variable)) + geom_point(shape=4, size=3) + geom_line() + theme_classic() + ylim(c(0,0.15)) + ylab('Methylation') + xlab('CpG position')
dev.off()
```
