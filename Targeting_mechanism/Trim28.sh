# Trim28 ChIP signal at DMRs
deeptools/bin/computeMatrix scale-regions \
-S WT_rep1_Trim28.bw WT_rep2_Trim28.bw TKOlike_rep1_Trim28.bw TKOlike_rep2_Trim28.bw DKOzero_rep1_Trim28.bw DKOzero_rep2_Trim28.bw \
-R ../DKO_TKO_DMRs_numerated_v3.bed \
-o DKO_TKO_DMRs_Trim28.tab.gz \
-a 15000 -b 15000 \
--binSize 50

deeptools/bin/plotHeatmap \
-m DKO_TKO_DMRs_Trim28.tab.gz \
-out figures/DKO_TKO_DMRs_Trim28_v3_profile.pdf \
--startLabel DMR \
--endLabel DMR \
--zMin 0 \
--zMax 0.25 \
--yMin -0.1 \
--yMax 1.5 \
--missingDataColor 'white' \
--colorMap 'Greys' \
--heatmapHeight 14
