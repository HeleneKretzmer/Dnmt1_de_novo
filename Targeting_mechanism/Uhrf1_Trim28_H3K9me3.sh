# Complex heatmap DMRs

```R
require(data.table)
require(ggplot2)
require(reshape2)
require(RColorBrewer)
require(EnrichedHeatmap)
require(rtracklayer)
require(circlize)

# Load data
TKO_track <- import('TKOlike.bw')
DKO_track <- import('DKOzero_P15.bw')
ERVK_track <- import('ERVK.bed')

Uhrf1_track <- import('Uhrf1_FE.bw')
Trim28_track <- import('Trim28_FE.bw')
H3K9me3_track <- import('H3K9me3_FE.bw')


# DMRs
DMRs <- import('../DKO_TKO_DMRs_numerated_v3.bed')

TKO_DMRs <- normalizeToMatrix(TKO_track, DMRs, value_column = 'score', extend = c(5000,5000), mean_mode = 'absolute', w = 250, background = NA, target_ratio = 0.25)
DKO_DMRs <- normalizeToMatrix(DKO_track, DMRs, value_column = 'score', extend = c(5000,5000), mean_mode = 'absolute', w = 250, background = NA, target_ratio = 0.25)
Uhrf1_DMRs <- normalizeToMatrix(Uhrf1_track, DMRs, value_column = 'score', extend = c(5000,5000), mean_mode = 'coverage', w = 250, background = 0, target_ratio = 0.25)
Trim28_DMRs <- normalizeToMatrix(Trim28_track, DMRs, value_column = 'score', extend = c(5000,5000), mean_mode = 'coverage', w = 250, background = 0, target_ratio = 0.25)
H3K9me3_DMRs <- normalizeToMatrix(H3K9me3_track, DMRs, value_column = 'score', extend = c(5000,5000), mean_mode = 'coverage', w = 250, background = 0, target_ratio = 0.25)
ERVK_DMRs <- normalizeToMatrix(ERVK_track, DMRs, value_column = 'name', extend = c(5000,5000), target_ratio = 0.25)

# heatmap
col_fun_meth <- colorRamp2(c(0, 0.15, 0.3, 0.5), c('blue', 'white', 'red', 'red'))
col_fun_Uhrf1 <- colorRamp2(c(0, 2), c('white', 'royalblue'))
col_fun_Trim28 <- colorRamp2(c(0, 10), c('white', 'darkolivegreen3'))
col_fun_H3K9me3 <- colorRamp2(c(0, 10), c('white', 'tomato'))
col_ERVK = c(
    "IAPLTR1a_Mm" = "goldenrod2",
    "IAPEY2_LTR"  = "gray76",
    "IAPEy-int"   = "gray60",
    "IAPEY3_LTR"  = "gray68",
    "IAP-d-int"   = "gray80",
    "IAPLTR2b"    = "gray52",
    "IAPLTR2_Mm"  = "goldenrod3",
    "IAPEz-int"   = "gold",
    "IAPEY3-int"  = "gray72",
    "IAPLTR4"     = "gray40",
    "IAPLTR3-int" = "gray44",
    "IAPLTR1_Mm"  = "goldenrod1",
    "IAPLTR4_I"   = "gray36",
    "IAPLTR3"     = "gray48",
    "IAPLTR2a"    = "gray56",
    "IAPEY_LTR"   = "gray64"
)

pdf('figures/ChIP_enrichment_DMRs_DKO.pdf', width = 20)
EnrichedHeatmap(DKO_DMRs, use_raster=T, col = col_fun_meth, axis_name_rot=90, column_title = 'DKO', name='DKO me', show_heatmap_legend=T, top_annotation = HeatmapAnnotation(enrich=anno_enriched(show_heatmap_legend=F, gp=gpar(col='black', lwd=3), ylim=c(0,0.3))), width=unit(4,'cm')) +
EnrichedHeatmap(Uhrf1_DMRs, use_raster=T, col = col_fun_Uhrf1, axis_name_rot=90, column_title = 'Uhrf1', name='Uhrf1 signal', show_heatmap_legend=T, top_annotation = HeatmapAnnotation(enrich=anno_enriched(show_heatmap_legend=F, gp=gpar(col='black', lwd=3), ylim=c(0.5,1.5))), width=unit(4,'cm')) +
EnrichedHeatmap(Trim28_DMRs, use_raster=T, col = col_fun_Trim28, axis_name_rot=90, column_title = 'Trim28', name='Trim28 signal', show_heatmap_legend=T, top_annotation = HeatmapAnnotation(enrich=anno_enriched(show_heatmap_legend=F, gp=gpar(col='black', lwd=3), ylim=c(0,6))), width=unit(4,'cm')) +
EnrichedHeatmap(H3K9me3_DMRs, use_raster=T, col = col_fun_H3K9me3, axis_name_rot=90, column_title = 'H3K9me3', name='H3K9me3 signal', show_heatmap_legend=T, top_annotation = HeatmapAnnotation(enrich=anno_enriched(show_heatmap_legend=F, gp=gpar(col='black', lwd=3), ylim=c(0,14))), width=unit(4,'cm')) +
EnrichedHeatmap(ERVK_DMRs, use_raster=T, col = col_ERVK, axis_name_rot=90, column_title = 'ERVK', name='ERVK', show_heatmap_legend=T, top_annotation = HeatmapAnnotation(enrich=anno_enriched(show_heatmap_legend=F, gp=gpar(col='black', lwd=3), ylim=c(0,0.6))), width=unit(4,'cm'))
dev.off()

```
