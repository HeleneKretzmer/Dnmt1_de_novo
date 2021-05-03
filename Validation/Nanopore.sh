# Unique vs multiple alignements
source virtual/nanopype/bin/activate

## mESC
samtools view -b -F 4 -@ 15 WT.bam | bedtools bamtobed -i stdin | python3 nanopype/rules/utils/alignment_cov.py mm9.chrom.sizes >WT.cov.bedGraph
samtools view -b -F 4 -F 256 -@ 15 WT.bam | bedtools bamtobed -i stdin | python3 nanopype/rules/utils/alignment_cov.py mm9.chrom.sizes >WT.cov.unique.bedGraph
## TKO
samtools view -b -F 4 -@ 15 TKO.bam | bedtools bamtobed -i stdin | python3 nanopype/rules/utils/alignment_cov.py mm9.chrom.sizes >TKO.cov.bedGraph
samtools view -b -F 4 -F 256 -@ 15 TKO.bam | bedtools bamtobed -i stdin | python3 nanopype/rules/utils/alignment_cov.py mm9.chrom.sizes >TKO.cov.unique.bedGraph
## DKO P15
samtools view -b -F 4 -@ 15 DKO.bam | bedtools bamtobed -i stdin | python3 nanopype/rules/utils/alignment_cov.py mm9.chrom.sizes >DKO.cov.bedGraph
samtools view -b -F 4 -F 256 -@ 15 DKO.bam | bedtools bamtobed -i stdin | python3 nanopype/rules/utils/alignment_cov.py mm9.chrom.sizes >DKO.cov.unique.bedGraph

## bigWig
for i in *.bedGraph; do
  o=`basename $i | sed 's/.bedGraph/.bw/'`
  echo $o;
  bedGraphToBigWig $i mm9.chrom.sizes $o;
done

## mean CpG coverage per DMR normalized to mean matching CR
for i in WT.cov DKO.cov TKO.cov; do
  bigWigAverageOverBed ${i}.bw ../DKO_TKO_DMRs_numerated_v1.bed stdout | cut -f1,5 >${i}.DMRs.tsv
  bigWigAverageOverBed ${i}.bw ../DKO_TKO_CRs_numerated_v1.bed stdout | sed 's/;iteration_[0-9]*//' | sort | bedtools groupby -g 1 -c 5 -o mean >${i}.CRs.tsv
  bigWigAverageOverBed ${i}.unique.bw ../DKO_TKO_DMRs_numerated_v1.bed stdout | cut -f1,5 >${i}.unique.DMRs.tsv
  bigWigAverageOverBed ${i}.unique.bw ../DKO_TKO_CRs_numerated_v1.bed stdout | sed 's/;iteration_[0-9]*//' | sort | bedtools groupby -g 1 -c 5 -o mean >${i}.unique.CRs.tsv
done

for s in DKO.hac.mm9.cov.bw  TKO.hac.mm9.cov.bw	WT.hac.mm9.cov.bw; do
  i=`basename ${s} | sed 's/.hac.mm9.cov.bw/.cov/'`;
  echo ${i};
  bigWigAverageOverBed ${s} ../DKO_TKO_DMRs_numerated_v1.bed stdout | cut -f1,5 >Nanopore.${i}.DMRs.tsv
  bigWigAverageOverBed ${s} ../DKO_TKO_CRs_numerated_v1.bed stdout | sed 's/;iteration_[0-9]*//' | sort | bedtools groupby -g 1 -c 5 -o mean >Nanopore.${i}.CRs.tsv
done

```R
require(ggplot2)

WT.cov.unique.CRs <- read.table('WT.cov.unique.CRs.tsv', header=F, col.names=c('region','CR'))
WT.cov.unique.DMRs <- read.table('WT.cov.unique.DMRs.tsv', header=F, col.names=c('region','DMR'))
WT.cov.unique <- merge(WT.cov.unique.DMRs, WT.cov.unique.CRs, by='region')
WT.cov.unique$FC <- WT.cov.unique$DMR/WT.cov.unique$CR
WT.cov.CRs <- read.table('WT.cov.CRs.tsv', header=F, col.names=c('region','CR'))
WT.cov.DMRs <- read.table('WT.cov.DMRs.tsv', header=F, col.names=c('region','DMR'))
WT.cov <- merge(WT.cov.DMRs, WT.cov.CRs, by='region')
WT.cov$FC <- WT.cov$DMR/WT.cov$CR

DKO.cov.unique.CRs <- read.table('DKO.cov.unique.CRs.tsv', header=F, col.names=c('region','CR'))
DKO.cov.unique.DMRs <- read.table('DKO.cov.unique.DMRs.tsv', header=F, col.names=c('region','DMR'))
DKO.cov.unique <- merge(DKO.cov.unique.DMRs, DKO.cov.unique.CRs, by='region')
DKO.cov.unique$FC <- DKO.cov.unique$DMR/DKO.cov.unique$CR
DKO.cov.CRs <- read.table('DKO.cov.CRs.tsv', header=F, col.names=c('region','CR'))
DKO.cov.DMRs <- read.table('DKO.cov.DMRs.tsv', header=F, col.names=c('region','DMR'))
DKO.cov <- merge(DKO.cov.DMRs, DKO.cov.CRs, by='region')
DKO.cov$FC <- DKO.cov$DMR/DKO.cov$CR

TKO.cov.unique.CRs <- read.table('TKO.cov.unique.CRs.tsv', header=F, col.names=c('region','CR'))
TKO.cov.unique.DMRs <- read.table('TKO.cov.unique.DMRs.tsv', header=F, col.names=c('region','DMR'))
TKO.cov.unique <- merge(TKO.cov.unique.DMRs, TKO.cov.unique.CRs, by='region')
TKO.cov.unique$FC <- TKO.cov.unique$DMR/TKO.cov.unique$CR
TKO.cov.CRs <- read.table('TKO.cov.CRs.tsv', header=F, col.names=c('region','CR'))
TKO.cov.DMRs <- read.table('TKO.cov.DMRs.tsv', header=F, col.names=c('region','DMR'))
TKO.cov <- merge(TKO.cov.DMRs, TKO.cov.CRs, by='region')
TKO.cov$FC <- TKO.cov$DMR/TKO.cov$CR

Nanopore.WT.cov.CRs <- read.table('Nanopore.WT.cov.CRs.tsv', header=F, col.names=c('region','CR'))
Nanopore.WT.cov.DMRs <- read.table('Nanopore.WT.cov.DMRs.tsv', header=F, col.names=c('region','DMR'))
Nanopore.WT.cov <- merge(Nanopore.WT.cov.DMRs, Nanopore.WT.cov.CRs, by='region')
Nanopore.WT.cov$FC <- Nanopore.WT.cov$DMR/Nanopore.WT.cov$CR

Nanopore.DKO.cov.CRs <- read.table('Nanopore.DKO.cov.CRs.tsv', header=F, col.names=c('region','CR'))
Nanopore.DKO.cov.DMRs <- read.table('Nanopore.DKO.cov.DMRs.tsv', header=F, col.names=c('region','DMR'))
Nanopore.DKO.cov <- merge(Nanopore.DKO.cov.DMRs, Nanopore.DKO.cov.CRs, by='region')
Nanopore.DKO.cov$FC <- Nanopore.DKO.cov$DMR/Nanopore.DKO.cov$CR

Nanopore.TKO.cov.CRs <- read.table('Nanopore.TKO.cov.CRs.tsv', header=F, col.names=c('region','CR'))
Nanopore.TKO.cov.DMRs <- read.table('Nanopore.TKO.cov.DMRs.tsv', header=F, col.names=c('region','DMR'))
Nanopore.TKO.cov <- merge(Nanopore.TKO.cov.DMRs, Nanopore.TKO.cov.CRs, by='region')
Nanopore.TKO.cov$FC <- Nanopore.TKO.cov$DMR/Nanopore.TKO.cov$CR


Nanopore.WT.cov$data <- 'Nanopore_WT'
Nanopore.DKO.cov$data <- 'Nanopore_DKO'
Nanopore.TKO.cov$data <- 'Nanopore_TKO'
WT.cov$data <- 'WT'
DKO.cov$data <- 'DKO'
TKO.cov$data <- 'TKO'
WT.cov.unique$data <- 'WT.unique'
DKO.cov.unique$data <- 'DKO.unique'
TKO.cov.unique$data <- 'TKO.unique'

data <- rbind(Nanopore.WT.cov,Nanopore.DKO.cov,Nanopore.TKO.cov, WT.cov,DKO.cov,TKO.cov, WT.cov.unique,DKO.cov.unique,TKO.cov.unique)

pdf('figures/FC_DMR_CR_coverage.pdf')
ggplot(data, aes(x=data, y=log2(FC))) + geom_boxplot(outlier.shape=NA) + theme_classic() + geom_hline(yintercept=0, color='grey')
dev.off()
```

# Exclude DMRs based on Nanopore sequencing
bedtools intersect -wa -wb -a DKO_TKO_DMRs_numerated_v2.bed -b DKO_TKO_diff.mm9.bedGraph | cut -f1,2,3,4,8 | bedtools groupby -g 1,2,3,4 -c 5 -o mean,median,max,count >DMRs_v2_diff_Nanopore.bed
bedtools intersect -wa -wb -a DKO_TKO_DMRs_numerated_v2.bed -b DKOzero_P15_TKO_diff.bed | cut -f1,2,3,4,8 | bedtools groupby -g 1,2,3,4 -c 5 -o mean,median,max,count >DMRs_v2_diff_Illumina.bed

```R
require(ggplot2)

Illumina <- read.table('DMRs_v2_diff_Illumina.bed')
Nanopore <- read.table('DMRs_v2_diff_Nanopore.bed')

colnames(Illumina) <- c('chr','start','end','ID','I_mean','I_median','I_max','I_count')
colnames(Nanopore) <- c('chr','start','end','ID','M_mean','M_median','M_max','M_count')

data <- merge(Illumina, Nanopore, by=c('chr','start','end','ID'))
data$M_diff_cat <- ifelse(data$M_mean > 0.05, 'strong', ifelse(data$M_mean < 0, 'none', 'weak'))
data <- subset(data, I_count>=10 & M_count>=5)

pdf('figures/Nanopore_exclude_DMRs_by_diff.pdf')
ggplot(data, aes(x=I_mean, y=M_mean)) + geom_point(aes(color=M_diff_cat)) + theme_classic() + xlim(-0.5,0.5) + ylim(-0.5,0.5) + geom_hline(yintercept=0) + geom_vline(xintercept=0) + scale_color_manual(values=c('grey90','black','grey50'))
dev.off()
```

# Methylation profile at DMRs
deeptools/bin/computeMatrix scale-regions \
-S DKO.hac.15x.mm9.bw TKO.hac.15x.mm9.bw \
-R DKO_TKO_DMRs_numerated_v3.bed \
-o DKO_TKO_DMRs.tab.gz \
-a 5000 -b 5000 \
--binSize 50

deeptools/bin/plotProfile \
-m DKO_TKO_DMRs.tab.gz \
-out figures/DKO_TKO_DMRs_numerated_v3_profile.pdf \
--plotType std \
--yMin -0.08 \
--yMax 0.55 \
--startLabel DMR \
--endLabel DMR \
--perGroup \
--colors 'orange' 'green'


# Heatmap DMRs
bedtools intersect -wa -wb -a ../DKO_TKO_DMRs_numerated_v3.bed -b DKO.hac.15x.mm9.bedGraph | bedtools groupby -g 1,2,3,4 -c 8 -o count | perl -ane 'if($F[-1]>=10){print $_}' | cut -f1,2,3,4 >DMRs_10CpG_DKO.bed

```R
require(data.table)
require(reshape2)
require(RColorBrewer)
require(EnrichedHeatmap)
require(rtracklayer)
require(circlize)

# Load data
DKO_track <- import('DKO.hac.15x.mm9.bw')

# DMRs
DMRs <- import('DMRs_10CpG_DKO.bed')
DKO_DMRs <- normalizeToMatrix(DKO_track, DMRs, value_column = 'score', extend = c(5000,5000), mean_mode = 'absolute', w = 500, background = NA, target_ratio = 0.25)

# heatmap
col_fun_meth <- colorRamp2(c(0, 0.1, 0.2, 0.5), c('blue', 'white', 'red', 'red'))

pdf('figures/DMRs_heatmap.pdf', width = 20)
EnrichedHeatmap(DKO_DMRs, col = col_fun_meth, axis_name_rot= 90, column_title = 'DKO', name='DKO me', show_heatmap_legend=T, top_annotation = HeatmapAnnotation(enrich=anno_enriched(show_heatmap_legend=F, gp=gpar(col='black', lwd=3), ylim=c(0,0.3))), width=unit(8,'cm'))
dev.off()

```

# DMR + CR average methylation scatter DKO~TKO
bedtools intersect -wa -wb -a DKO_TKO_DMRs_numerated_v3.bed -b WT.hac.15x.mm9.bedGraph | cut -f1,2,3,4,8 | bedtools groupby -g 1,2,3,4 -c 5 -o mean,median,max,count >DMRs_mESCmeM.txt
bedtools intersect -wa -wb -a DKO_TKO_DMRs_numerated_v3.bed -b DKO.hac.15x.mm9.bedGraph | cut -f1,2,3,4,8 | bedtools groupby -g 1,2,3,4 -c 5 -o mean,median,max,count >DMRs_DKOmeM.txt
bedtools intersect -wa -wb -a DKO_TKO_DMRs_numerated_v3.bed -b TKO.hac.15x.mm9.bedGraph | cut -f1,2,3,4,8 | bedtools groupby -g 1,2,3,4 -c 5 -o mean,median,max,count >DMRs_TKOmeM.txt
bedtools intersect -wa -wb -a DKO_TKO_CRs_numerated_v3.bed -b WT.hac.15x.mm9.bedGraph | cut -f1,2,3,4,8 | bedtools groupby -g 1,2,3,4 -c 5 -o mean,median,max,count >CRs_mESCmeM.txt
bedtools intersect -wa -wb -a DKO_TKO_CRs_numerated_v3.bed -b DKO.hac.15x.mm9.bedGraph | cut -f1,2,3,4,8 | bedtools groupby -g 1,2,3,4 -c 5 -o mean,median,max,count >CRs_DKOmeM.txt
bedtools intersect -wa -wb -a DKO_TKO_CRs_numerated_v3.bed -b TKO.hac.15x.mm9.bedGraph | cut -f1,2,3,4,8 | bedtools groupby -g 1,2,3,4 -c 5 -o mean,median,max,count >CRs_TKOmeM.txt

```R
require(ggplot2)
require(vioplot)
require(plyr)

DMR_mESC <-  read.table('DMRs_mESCmeM.txt')
DMR_DKO <-  read.table('DMRs_DKOmeM.txt')
DMR_TKO <- read.table('DMRs_TKOmeM.txt')
CR_mESC <-  read.table('CRs_mESCmeM.txt')
CR_DKO <-  read.table('CRs_DKOmeM.txt')
CR_TKO <- read.table('CRs_TKOmeM.txt')

DMR <- Reduce(function(x, y) merge(x, y, all.x=T, all.y=T, by=c('V4')), list(DMR_mESC[,c(4,5)], DMR_DKO[,c(4,5)], DMR_TKO[,c(4,5)]))
colnames(DMR) <- c('ID','DMR_mESC','DMR_TKO','DMR_DKOzero')
CR <- Reduce(function(x, y) merge(x, y, all.x=T, all.y=T, by=c('V4')), list(CR_mESC[,c(4,5)], CR_DKO[,c(4,5)], CR_TKO[,c(4,5)]))
colnames(CR) <- c('ID','CR_mESC','CR_TKO','CR_DKOzero')

pdf('figures/DKO_TKO_DMRs_vioplot_avg_per_region.pdf', width=14)
vioplot(DMR_mESC$V5, DMR_DKO$V5, DMR_TKO$V5, CR_mESC$V5, CR_DKO$V5, CR_TKO$V5, names=c('DMR_mESC','DMR_DKO','DMR_TKO', 'CR_mESC','CR_DKO','CR_TKO'), col=c('grey','orange','forestgreen','grey','orange','forestgreen'))
dev.off()
```

# DMR and CR scatter plots
```R
require(ggplot2)
require(vioplot)
require(plyr)

DMR_DKO <-  read.table('DMRs_DKOmeM.txt')
DMR_TKO <- read.table('DMRs_TKOmeM.txt')
CR_DKO <-  read.table('CRs_DKOmeM.txt')
CR_TKO <- read.table('CRs_TKOmeM.txt')

DMRs <- merge(DMR_TKO[,4:5], DMR_DKO[,4:5], by='V4')
colnames(DMRs) <- c('ID','TKO','DKO')
CRs <- merge(CR_TKO[,4:5], CR_DKO[,4:5], by='V4')
colnames(CRs) <- c('ID','TKO','DKO')

set.seed(12)
subCRs <- CRs[sample(nrow(CRs), nrow(DMRs)), ]


pdf('figures/DMRs_CRs_scatter.pdf')
ggplot(DMRs, aes(x = TKO, y=DKO)) + geom_point(aes(color=ifelse(DMRs$DKO - DMRs$TKO < 0.35, DMRs$DKO - DMRs$TKO, 0.35))) + theme_classic() + xlim(c(0,1)) + ylim(c(0,1)) + scale_color_gradient2(limits=c(-0.35,0.35), midpoint=0, low='darkblue', mid='gray90', high='firebrick3', space='Lab') + xlab('Mean methylation TKO') + ylab('Mean methylation DKO') + ggtitle('DMRs') + geom_abline(slope=1, color='grey80') + labs(color = 'DKO-TKO')
ggplot(subCRs, aes(x = TKO, y=DKO)) + geom_point(aes(color=ifelse(subCRs$DKO - subCRs$TKO < 0.35, subCRs$DKO - subCRs$TKO, 0.35))) + theme_classic() + xlim(c(0,1)) + ylim(c(0,1)) + scale_color_gradient2(limits=c(-0.35,0.35), midpoint=0, low='darkblue', mid='gray90', high='firebrick3', space='Lab') + xlab('Mean methylation TKO') + ylab('Mean methylation DKO') + ggtitle('CRs (sampled)') + geom_abline(slope=1, color='grey80') + labs(color = 'DKO-TKO')

ggplot(DMRs, aes(x = TKO, y=DKO)) + geom_point(aes(color=ifelse(DMRs$DKO - DMRs$TKO < 0.35, DMRs$DKO - DMRs$TKO, 0.35))) + theme_classic() + xlim(c(0,1)) + ylim(c(0,1)) + scale_color_gradient2(limits=c(-0.35,0.35), midpoint=0, low='darkblue', mid='gray90', high='firebrick3', space='Lab') + xlab('Mean methylation TKO') + ylab('Mean methylation DKO') + ggtitle('DMRs') + geom_abline(slope=1, color='grey80') + labs(color = 'DKO-TKO') + theme(legend.position='none')
ggplot(subCRs, aes(x = subCRs$TKO, y=subCRs$DKO)) + geom_point(aes(color=ifelse(subCRs$DKO - subCRs$TKO < 0.35, subCRs$DKO - subCRs$TKO, 0.35))) + theme_classic() + xlim(c(0,1)) + ylim(c(0,1)) + scale_color_gradient2(limits=c(-0.35,0.35), midpoint=0, low='darkblue', mid='gray90', high='firebrick3', space='Lab') + xlab('Mean methylation TKO') + ylab('Mean methylation DKO') + ggtitle('CRs (sampled)') + geom_abline(slope=1, color='grey80') + labs(color = 'DKO-TKO') + theme(legend.position='none')
dev.off()
```

# Avg. Illumina methylation in deleted vs present DMRs
bedtools intersect -wa -wb -a DKO_TKO_DMRs_numerated_v3.bed -b WT.bed | cut -f1,2,3,4,8 | bedtools groupby -g 1,2,3,4 -c 5 -o mean,median,max,count >DMRs_mESCme.txt
bedtools intersect -wa -wb -a DKO_TKO_DMRs_numerated_v3.bed -b DKOzero_P15.bed | cut -f1,2,3,4,8 | bedtools groupby -g 1,2,3,4 -c 5 -o mean,median,max,count >DMRs_DKOme.txt
bedtools intersect -wa -wb -a DKO_TKO_DMRs_numerated_v3.bed -b TKOlike.bed | cut -f1,2,3,4,8 | bedtools groupby -g 1,2,3,4 -c 5 -o mean,median,max,count >DMRs_TKOme.txt

bedtools intersect -wa -wb -a DKO_TKO_DMRs_deleted.bed -b WT.bed | cut -f1,2,3,4,8 | bedtools groupby -g 1,2,3,4 -c 5 -o mean,median,max,count >Deletions_mESCme.txt
bedtools intersect -wa -wb -a DKO_TKO_DMRs_deleted.bed -b DKOzero_P15.bed | cut -f1,2,3,4,8 | bedtools groupby -g 1,2,3,4 -c 5 -o mean,median,max,count >Deletions_DKOme.txt
bedtools intersect -wa -wb -a DKO_TKO_DMRs_deleted.bed -b TKOlike.bed | cut -f1,2,3,4,8 | bedtools groupby -g 1,2,3,4 -c 5 -o mean,median,max,count >Deletions_TKOme.txt

```R
require(ggplot2)
require(vioplot)
require(plyr)

DMR_mESC <-  read.table('DMRs_mESCme.txt')
DMR_DKO <-  read.table('DMRs_DKOme.txt')
DMR_TKO <- read.table('DMRs_TKOme.txt')
Del_mESC <-  read.table('Deletions_mESCme.txt')
Del_DKO <-  read.table('Deletions_DKOme.txt')
Del_TKO <- read.table('Deletions_TKOme.txt')

remaining <- Reduce(function(x, y) merge(x, y, all.x=T, all.y=T, by=c('V4')), list(DMR_mESC[,c(4,5)], DMR_TKO[,c(4,5)], DMR_DKO[,c(4,5)]))
colnames(remaining) <- c('ID','mESC','TKO','DKO')
exclude <- Reduce(function(x, y) merge(x, y, all.x=T, all.y=T, by=c('V4')), list(Del_mESC[,c(4,5)], Del_TKO[,c(4,5)], Del_DKO[,c(4,5)]))
colnames(exclude) <- c('ID','mESC','TKO','DKO')

pdf('figures/Illumina_DKO_TKO_DMRs_Deletions_vioplots_avg_per_region.pdf', width=10)
vioplot(DMR_mESC$V5, DMR_DKO$V5, DMR_TKO$V5, Del_mESC$V5, Del_DKO$V5, Del_TKO$V5, names=c('DMR_mESC','DMR_DKO','DMR_TKO', 'Del_mESC','Del_DKO','Del_TKO'), col=c('grey','orange','forestgreen','grey','orange','forestgreen'))
dev.off()
```

# Read-level analysis
```
require(ggplot2)
require(ggalluvial)

WT_DMRs <- read.table('WT_DMR_mean_methylation.bed', header=F, col.names=c('chr','start','end','DMR','ID','count','mean_methylation'))
J1_TKOlike_DMRs <- read.table('TKO_DMR_mean_methylation.bed', header=F, col.names=c('chr','start','end','DMR','ID','count','mean_methylation'))
J1_DKO_DMRs <- read.table('DKO_DMR_mean_methylation.bed', header=F, col.names=c('chr','start','end','DMR','ID','count','mean_methylation'))

WT_DMRs <- data.frame(table(cut(subset(WT_DMRs, count>=10)$mean_methylation, breaks=c(0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 1), include.lowest=T)))
J1_TKOlike_DMRs <- data.frame(table(cut(subset(J1_TKOlike_DMRs, count>=10)$mean_methylation, breaks=c(0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 1), include.lowest=T)))
J1_DKO_DMRs <- data.frame(table(cut(subset(J1_DKO_DMRs, count>=10)$mean_methylation, breaks=c(0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 1), include.lowest=T)))

WT_DMRs$frac <- 100*WT_DMRs$Freq/sum(WT_DMRs$Freq)
WT_DMRs$type <- 'WT'
J1_TKOlike_DMRs$frac <- 100*J1_TKOlike_DMRs$Freq/sum(J1_TKOlike_DMRs$Freq)
J1_TKOlike_DMRs$type <- 'TKO'
J1_DKO_DMRs$frac <- 100*J1_DKO_DMRs$Freq/sum(J1_DKO_DMRs$Freq)
J1_DKO_DMRs$type <- 'DKO'

data <- rbind(WT_DMRs, J1_TKOlike_DMRs, J1_DKO_DMRs)
data$type <- factor(data$type, levels=c('TKO', 'DKO', 'WT'))

pdf('figures/read_methylation_level_Nanopore.pdf')
ggplot(data, aes(x=type, y=frac, fill=Var1)) + geom_bar(stat='identity') + theme_classic() + scale_fill_manual(values=c('dodgerblue3','cadetblue3','chartreuse3','greenyellow','yellow','goldenrod1','orange','firebrick'))
dev.off()
```
