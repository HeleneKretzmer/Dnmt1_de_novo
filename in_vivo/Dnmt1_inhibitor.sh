# Normal heatmap IAPEz-int
for i in WT_E3.5.bed DKO_0.35uM_1i_E35.bed Dnmt1_KO_E3.5.bed WT_E6.5_Epi.bed DKO_0.35uM_1i_Epi.bed Dnmt1_KO_E6.5_Epi.bed; do
  n=`basename $i | sed 's/.bed//'`;
  echo ${n};
  bedtools intersect -wa -wb -a IAPEz-int.bed -b ${i} | bedtools groupby -g 1,2,3,4 -c 8 -o mean,count | perl -ane 'if($F[-1]>=5){print $_}' | cut -f1,2,3,5 >${n}_IAPEzint.tsv;
done

```R
require(ggplot2)
require(RColorBrewer)
require(ComplexHeatmap)
require(circlize)

# load data
WT_E35 <- read.table('/project/GhostPeaks/in_vivo/WT_E3.5_IAPEzint.tsv', header=F, col.names=c('chr','start','end','WT_E35'))
Dnmt3ab_D1i_E35 <- read.table('/project/GhostPeaks/in_vivo/DKO_0.35uM_1i_E35_IAPEzint.tsv', header=F, col.names=c('chr','start','end','Dnmt3ab_D1i_E35'))
Dnmt1_E35 <- read.table('/project/GhostPeaks/in_vivo/Dnmt1_KO_E3.5_IAPEzint.tsv', header=F, col.names=c('chr','start','end','Dnmt1_E35'))
WT_E65 <- read.table('/project/GhostPeaks/in_vivo/WT_E6.5_Epi_IAPEzint.tsv', header=F, col.names=c('chr','start','end','WT_E65'))
Dnmt3ab_D1i_E65 <- read.table('/project/GhostPeaks/in_vivo/DKO_0.35uM_1i_Epi_IAPEzint.tsv', header=F, col.names=c('chr','start','end','Dnmt3ab_D1i_E65'))
Dnmt1_E65 <- read.table('/project/GhostPeaks/in_vivo/Dnmt1_KO_E6.5_Epi_IAPEzint.tsv', header=F, col.names=c('chr','start','end','Dnmt1_E65'))

data <- Reduce(function(df1, df2) merge(df1, df2, by = c('chr','start','end'), all.x = TRUE), list(WT_E35,Dnmt3ab_D1i_E35,Dnmt1_E35, WT_E65,Dnmt3ab_D1i_E65,Dnmt1_E65))

pdf('figures/Heatmap_IAPEzint_E3.5_E6.5.pdf')
Heatmap(data[,c(4,9, 5,10, 6,11, 7,12, 8)], cluster_columns=F, show_row_names=F, col=colorRamp2(c(0, 0.15, 0.8, 1), c('blue', 'white', 'red', 'red')))
dev.off()
```

# Delta E3.5 vs E6.5 at DMR
deeptools/bin/computeMatrix scale-regions \
-S delta_WT.bw delta_Dnmt1.bw delta_Dnmt3abD1i.bw \
-R IAPEz-int.bed \
-o IAPEz-int_delta.tab.gz \
-a 10000 -b 10000 \
--binSize 100

deeptools/bin/plotProfile \
-m IAPEz-int_delta.tab.gz \
-out figures/IAPEz-int_delta_profile.pdf \
--startLabel IAPEz-int \
--endLabel IAPEz-int \
--yMin -0.75 \
--yMax 0.75 \
--plotType=std \
--perGroup


# DKO vs DKO Dnmt1i E35
deeptools/bin/computeMatrix scale-regions \
-S Dnmt3ab_DKO_E3.5.bw DKO_0.35uM_1i_E35.bw \
-R IAPEz-int.bed \
-o IAPEz-int_E35.tab.gz \
-a 10000 -b 10000 \
--binSize 100

deeptools/bin/plotProfile \
-m IAPEz-int_E35.tab.gz \
-out figures/IAPEz-int_E35_profile.pdf \
--startLabel IAPEz-int \
--endLabel IAPEz-int \
--yMin -0.75 \
--yMax 0.75 \
--plotType=std \
--perGroup


# DKO vs DKO Dnmt1i E65
deeptools/bin/computeMatrix scale-regions \
-S Dnmt3ab_DKO_Epi.bw DKO_0.35uM_1i_Epi.bw \
-R IAPEz-int.bed \
-o IAPEz-int_E65.tab.gz \
-a 10000 -b 10000 \
--binSize 100

deeptools/bin/plotProfile \
-m IAPEz-int_E35.tab.gz \
-out figures/IAPEz-int_E65_profile.pdf \
--startLabel IAPEz-int \
--endLabel IAPEz-int \
--yMin -0.75 \
--yMax 0.75 \
--plotType=std \
--perGroup
