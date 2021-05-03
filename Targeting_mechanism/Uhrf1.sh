# Uhrf1 KO methylation at DMRs
deeptools/bin/computeMatrix scale-regions \
-S DKOzero_P15.bw DKO_UHRF1_A10_P15.bw DKO_UHRF1_C5_P15.bw \
-R ../DKO_TKO_DMRs_numerated_v3.bed \
-o DKO_TKO_DMRs_v3.tab.gz \
-a 2500 -b 2500 \
--binSize 50

/project/bioinf_meissner/src/deeptools/bin/plotHeatmap \
-m DKO_TKO_DMRs_v3.tab.gz \
-out figures/DKO_TKO_DMRs_v3_profile.pdf \
--startLabel DMR \
--endLabel DMR \
--zMin 0 \
--zMax 1 \
--yMin -0.1 \
--yMax 0.6 \
--missingDataColor 'white' \
--colorMap 'Greys' \
--heatmapHeight 14

# Uhrf1 ChIP enrichment at DMRs
deeptools/bin/computeMatrix scale-regions \
-S WT_Uhrf1_rep1.bw WT_Uhrf1_rep2.bw TKOlike_Uhrf1_rep1.bw TKOlike_Uhrf1_rep2.bw DKOzero_Uhrf1_rep1.bw DKOzero_Uhrf1_rep2.bw \
-R ../DKO_TKO_DMRs_numerated_v3.bed \
-o DKO_TKO_DMRs_Uhrf1.tab.gz \
-a 15000 -b 15000 \
--binSize 50

deeptools/bin/plotHeatmap \
-m DKO_TKO_DMRs_Uhrf1.tab.gz \
-out figures/DKO_TKO_DMRs_Uhrf1_v3_profile.pdf \
--startLabel DMR \
--endLabel DMR \
--zMin 0 \
--zMax 0.25 \
--yMin -0.1 \
--yMax 0.25 \
--missingDataColor 'white' \
--colorMap 'Greys' \
--heatmapHeight 14
