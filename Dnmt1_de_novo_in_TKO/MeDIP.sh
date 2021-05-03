# Percent of reads in DMR/CR
for i in MeDIP_*.bam; do
    n=`basename $i | sed 's/.bam//' | sed 's/MeDIP_//'`
    echo $n
    # total reads
    samtools view -c -F 256 -F 4 ${i} | awk -v n=$n '{print "All\t"n"\t"$0}' >>MeDIP_read_stats.txt
    # reads in DMRs
    samtools view -c -F 256 -F 4 -L /project/GhostPeaks/DKO_TKO_DMRs_numerated_v3.bed ${i} | awk -v n=$n '{print "DMR\t"n"\t"$0}' >>MeDIP_read_stats.txt
    # reads in CRs
    samtools view -c -F 256 -F 4 -L /project/GhostPeaks/DKO_TKO_CRs_numerated_v3.bed ${i} | awk -v n=$n '{print "CR\t"n"\t"$0}' >>MeDIP_read_stats.txt
done

```R
require(ggplot2)
require(plyr)

data <- read.table('MeDIP_read_stats.txt')
colnames(data) <- c('region','ident','count')

data$ID <- c('p1_rep1','p1_rep2','p1_rep3','p5_rep1','p5_rep2','p5_rep3','p10_rep1','p10_rep2','p10_rep3')
data$P <- gsub('_.*','',data$ID)
data$P <- factor(data$P, levels=c('p1','p5','p10'))

data <- subset(data, region != 'All')
data$norm_count <- data$count / rep(c(1,1000),9)

pdf('figures/MeDIP_read_stats.pdf')
ggplot(data, aes(x=region, y=norm_count)) + geom_boxplot() + geom_point() + theme_classic() + facet_wrap(~P)
dev.off()
```

# Heatmap enrichment of signal at DMRs/CRs
deeptools/bin/computeMatrix scale-regions \
-S MeDip_P1_counts.bw MeDip_P5_counts.bw MeDip_P10_counts.bw \
-R ../DKO_TKO_DMRs_numerated_v3.bed \
-o DKO_TKO_DMRs.tab.gz \
-a 5000 -b 5000 \
--binSize 100 \
--missingDataAsZero

deeptools/bin/plotHeatmap \
-m DKO_TKO_DMRs.tab.gz \
-out figures/DKO_TKO_DMRs_numerated_v3_profile.pdf \
--startLabel DMR \
--endLabel DMR \
--zMin 0 \
--zMax 8 \
--missingDataColor 'white' \
--colorMap 'Greys' \
--heatmapHeight 14
