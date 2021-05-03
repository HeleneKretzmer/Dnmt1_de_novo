# H3K9m3 signal at DMRs
/project/bioinf_meissner/src/deeptools/bin/computeMatrix scale-regions \
-S DKOzero_p5_rep1_H3K9me3.bw DKOzero_p5_rep2_H3K9me3.bw TKOlike_rep1_H3K9me3.bw TKOlike_rep2_H3K9me3.bw V65_rep1_H3K9me3.bw V65_rep2_H3K9me3.bw \
-R ../../../DKO_TKO_DMRs_numerated_v3.bed \
-o /project/GhostPeaks/Overview/DKO_TKO_DMRs_H3K9me3.tab.gz \
-a 15000 -b 15000 \
--binSize 50

/project/bioinf_meissner/src/deeptools/bin/plotHeatmap \
-m /project/GhostPeaks/Overview/DKO_TKO_DMRs_H3K9me3.tab.gz \
-out /project/GhostPeaks/Overview/figures/DKO_TKO_DMRs_H3K9me3_v3_profile.pdf \
--startLabel DMR \
--endLabel DMR \
--zMin 0 \
--zMax 1.5 \
--yMin -0.1 \
--yMax 2 \
--missingDataColor 'white' \
--colorMap 'Greys' \
--heatmapHeight 14


# Overlap DMRs with H3K9me3 peaks
bedtools intersect -u -a ../DKO_TKO_DMRs_numerated_v3.bed -b TKOlike_H3K9me3.broadPeak.bed | wc -l
1340
wc -l ../DKO_TKO_DMRs_numerated_v3.bed TKOlike_H3K9me3.broadPeak.bed
   1515 ../DKO_TKO_DMRs_numerated_v3.bed
 104618 TKOlike_H3K9me3.broadPeak.bed

bedtools intersect -u -a ../DKO_TKO_DMRs_numerated_v3.bed -b ../data/ChIP/peaks/DKOzero_H3K9me3.broadPeak.bed | wc -l
1359
wc -l ../DKO_TKO_DMRs_numerated_v3.bed DKOzero_H3K9me3.broadPeak.bed
   1515 ../DKO_TKO_DMRs_numerated_v3.bed
  92462 DKOzero_H3K9me3.broadPeak.bed

```R
library(VennDiagram)

pdf('figures/DKO_TKO_DMRs_v3_H3K9me3_venn.pdf')
grid.newpage()
draw.pairwise.venn(area1      = 104618,
                                area2      = 1515,
                                cross.area = 1340,
                                category   = c('TKO H3K9me3', 'DMRs'))

grid.newpage()
draw.pairwise.venn(area1      = 92462,
                                area2      = 1515,
                                cross.area = 1359,
                                category   = c('DKO H3K9me3', 'DMRs'))
dev.off()
```



# Vioplot methylation at H3K9me3 peaks
bedtools intersect -u -a DKOzero_P15.bed -b DKOzero_H3K9me3.broadPeak.bed >DKOzero_H3K9me3_DKOme.bed
bedtools intersect -v -a DKOzero_P15.bed -b DKOzero_H3K9me3.broadPeak.bed >DKOzero_notH3K9me3_DKOme.bed

```R
require(vioplot)

K9 <- read.table('DKOzero_H3K9me3_DKOme.bed', header=F)
notK9 <- read.table('DKOzero_notH3K9me3_DKOme.bed', header=F)

pdf('figures/DKOme_at_H3K9me3_peaks_genomewide_vioplots.pdf',width=10)
vioplot(list(K9[,4], notK9[,4]), names=c('H3K9me3 peaks','Not H3K9me3 peaks'), col=c('orange', 'orange'))
dev.off()
```
