grep -v nan WT_E35_readlevel.bed | bedtools intersect -wa -wb -a ../in_vivo/IAPEz-int.bed -b stdin >WT_E3.5_IAP_readlevel.bed
grep -v nan WT_E35_readlevel.bed | bedtools intersect -wa -wb -a ../DKO_TKO_CRs_numerated_v3.bed -b stdin >WT_E3.5_CR_readlevel.bed
grep -v nan Dnmt1_E35_readlevel.bed | bedtools intersect -wa -wb -a ../in_vivo/IAPEz-int.bed -b stdin >Dnmt1_KO_E3.5_IAP_readlevel.bed
grep -v nan Dnmt1_E35_readlevel.bed | bedtools intersect -wa -wb -a ../DKO_TKO_CRs_numerated_v3.bed -b stdin >Dnmt1_KO_E3.5_CR_readlevel.bed
grep -v nan DKOD1i_E35_readlevel.bed | bedtools intersect -wa -wb -a ../in_vivo/IAPEz-int.bed -b stdin >DKO_0.35uM_1i_E35_IAP_readlevel.bed
grep -v nan DKOD1i_E35_readlevel.bed | bedtools intersect -wa -wb -a ../DKO_TKO_CRs_numerated_v3.bed -b stdin >DKO_0.35uM_1i_E35_CR_readlevel.bed
grep -v nan WT_E65_readlevel.bed | bedtools intersect -wa -wb -a ../in_vivo/IAPEz-int.bed -b stdin >WT_E6.5_IAP_readlevel.bed
grep -v nan WT_E65_readlevel.bed | bedtools intersect -wa -wb -a ../DKO_TKO_CRs_numerated_v3.bed -b stdin >WT_E6.5_CR_readlevel.bed
grep -v nan Dnmt1_E65_readlevel.bed | bedtools intersect -wa -wb -a ../in_vivo/IAPEz-int.bed -b stdin >DKO_0.35uM_1i_Epi_IAP_readlevel.bed
grep -v nan Dnmt1_E65_readlevel.bed | bedtools intersect -wa -wb -a ../DKO_TKO_CRs_numerated_v3.bed -b stdin >DKO_0.35uM_1i_Epi_CR_readlevel.bed
grep -v nan DKOD1i_E65_readlevel.bed | bedtools intersect -wa -wb -a ../in_vivo/IAPEz-int.bed -b stdin >Dnmt1_KO_E6.5_Epi_IAP_readlevel.bed
grep -v nan DKOD1i_E65_readlevel.bed | bedtools intersect -wa -wb -a ../DKO_TKO_CRs_numerated_v3.bed -b stdin >Dnmt1_KO_E6.5_Epi_CR_readlevel.bed



```R
require(ggplot2)
require(ggalluvial)

WT_E3.5_DMRs <- read.table('WT_E3.5_IAP_readlevel.bed', header=F, col.names=c('chr','start','end','DMR','c','s','e','CpG_pattern','n_cpgs','n_cpgs_methyl','discordance_score','mean_methylation'))
Dnmt1_KO_E3.5_DMRs <- read.table('Dnmt1_KO_E3.5_IAP_readlevel.bed', header=F, col.names=c('chr','start','end','DMR','c','s','e','CpG_pattern','n_cpgs','n_cpgs_methyl','discordance_score','mean_methylation'))
DKO_0.35uM_1i_E35_DMRs <- read.table('DKO_0.35uM_1i_E35_IAP_readlevel.bed', header=F, col.names=c('chr','start','end','DMR','c','s','e','CpG_pattern','n_cpgs','n_cpgs_methyl','discordance_score','mean_methylation'))
WT_E6.5_Epi_DMRs <- read.table('WT_E6.5_IAP_readlevel.bed', header=F, col.names=c('chr','start','end','DMR','c','s','e','CpG_pattern','n_cpgs','n_cpgs_methyl','discordance_score','mean_methylation'))
Dnmt1_KO_E6.5_Epi_DMRs <- read.table('Dnmt1_KO_E6.5_Epi_IAP_readlevel.bed', header=F, col.names=c('chr','start','end','DMR','c','s','e','CpG_pattern','n_cpgs','n_cpgs_methyl','discordance_score','mean_methylation'))
DKO_0.35uM_1i_Epi_DMRs <- read.table('DKO_0.35uM_1i_Epi_IAP_readlevel.bed', header=F, col.names=c('chr','start','end','DMR','c','s','e','CpG_pattern','n_cpgs','n_cpgs_methyl','discordance_score','mean_methylation'))

WT_E3.5 <- data.frame(table(cut(WT_E3.5_DMRs$mean_methylation, breaks=c(0, 0.0000001, 0.25, 0.5, 0.75, 1), include.lowest=T)))
WT_E6.5_Epi <- data.frame(table(cut(WT_E6.5_Epi_DMRs$mean_methylation, breaks=c(0, 0.0000001, 0.25, 0.5, 0.75, 1), include.lowest=T)))
Dnmt1_KO_E3.5 <- data.frame(table(cut(Dnmt1_KO_E3.5_DMRs$mean_methylation, breaks=c(0, 0.0000001, 0.25, 0.5, 0.75, 1), include.lowest=T)))
Dnmt1_KO_E6.5_Epi <- data.frame(table(cut(Dnmt1_KO_E6.5_Epi_DMRs$mean_methylation, breaks=c(0, 0.0000001, 0.25, 0.5, 0.75, 1), include.lowest=T)))
DKO_0.35uM_1i_E35 <- data.frame(table(cut(DKO_0.35uM_1i_E35_DMRs$mean_methylation, breaks=c(0, 0.0000001, 0.25, 0.5, 0.75, 1), include.lowest=T)))
DKO_0.35uM_1i_Epi <- data.frame(table(cut(DKO_0.35uM_1i_Epi_DMRs$mean_methylation, breaks=c(0, 0.0000001, 0.25, 0.5, 0.75, 1), include.lowest=T)))

WT_E3.5$frac <- 100*WT_E3.5$Freq/sum(WT_E3.5$Freq)
WT_E3.5$type <- 'E3.5'
WT_E3.5$group <- 'WT'
WT_E6.5_Epi$frac <- 100*WT_E6.5_Epi$Freq/sum(WT_E6.5_Epi$Freq)
WT_E6.5_Epi$type <- 'E6.5_Epi'
WT_E6.5_Epi$group <- 'WT'
Dnmt1_KO_E3.5$frac <- 100*Dnmt1_KO_E3.5$Freq/sum(Dnmt1_KO_E3.5$Freq)
Dnmt1_KO_E3.5$type <- 'E3.5'
Dnmt1_KO_E3.5$group <- 'Dnmt1_KO'
Dnmt1_KO_E6.5_Epi$frac <- 100*Dnmt1_KO_E6.5_Epi$Freq/sum(Dnmt1_KO_E6.5_Epi$Freq)
Dnmt1_KO_E6.5_Epi$type <- 'E6.5_Epi'
Dnmt1_KO_E6.5_Epi$group <- 'Dnmt1_KO'
DKO_0.35uM_1i_E35$frac <- 100*DKO_0.35uM_1i_E35$Freq/sum(DKO_0.35uM_1i_E35$Freq)
DKO_0.35uM_1i_E35$type <- 'E3.5'
DKO_0.35uM_1i_E35$group <- 'DKO_0.35uM_1i'
DKO_0.35uM_1i_Epi$frac <- 100*DKO_0.35uM_1i_Epi$Freq/sum(DKO_0.35uM_1i_Epi$Freq)
DKO_0.35uM_1i_Epi$type <- 'E6.5_Epi'
DKO_0.35uM_1i_Epi$group <- 'DKO_0.35uM_1i'

data <- rbind(WT_E3.5, WT_E6.5_Epi, Dnmt1_KO_E3.5, Dnmt1_KO_E6.5_Epi, DKO_0.35uM_1i_E35, DKO_0.35uM_1i_Epi)
data$type <- factor(data$type, levels=c('E3.5', 'E6.5_Epi'))
data$Var1 <- factor(data$Var1, levels=rev(levels(data$Var1)))

pdf('read_methylation_level_in_vivo.pdf')
ggplot(data, aes(x = type, stratum = Var1, alluvium = Var1, y = frac, fill = Var1, label = Var1)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow(alpha = .5) +
  geom_stratum(color = 'lightgrey') +
  theme_classic() +
  scale_fill_manual(values=rev(c('dodgerblue3','chartreuse3','yellow','orange','firebrick'))) + 
  facet_grid(~group)
dev.off()



Dnmt1_KO_E3.5_CRs <- read.table('Dnmt1_KO_E3.5_CR_readlevel.bed', header=F, col.names=c('chr','start','end','DMR','c','s','e','CpG_pattern','n_cpgs','n_cpgs_methyl','discordance_score','mean_methylation'))
DKO_0.35uM_1i_E35_CRs <- read.table('DKO_0.35uM_1i_E35_CR_readlevel.bed', header=F, col.names=c('chr','start','end','DMR','c','s','e','CpG_pattern','n_cpgs','n_cpgs_methyl','discordance_score','mean_methylation'))
Dnmt1_KO_E6.5_Epi_CRs <- read.table('Dnmt1_KO_E6.5_Epi_CR_readlevel.bed', header=F, col.names=c('chr','start','end','DMR','c','s','e','CpG_pattern','n_cpgs','n_cpgs_methyl','discordance_score','mean_methylation'))
DKO_0.35uM_1i_Epi_CRs <- read.table('DKO_0.35uM_1i_Epi_CR_readlevel.bed', header=F, col.names=c('chr','start','end','DMR','c','s','e','CpG_pattern','n_cpgs','n_cpgs_methyl','discordance_score','mean_methylation'))

Dnmt1_KO_E3.5 <- data.frame(table(cut(Dnmt1_KO_E3.5_CRs$mean_methylation, breaks=c(0, 0.0000001, 0.25, 0.5, 0.75, 1), include.lowest=T)))
Dnmt1_KO_E6.5_Epi <- data.frame(table(cut(Dnmt1_KO_E6.5_Epi_CRs$mean_methylation, breaks=c(0, 0.0000001, 0.25, 0.5, 0.75, 1), include.lowest=T)))
DKO_0.35uM_1i_E35 <- data.frame(table(cut(DKO_0.35uM_1i_E35_CRs$mean_methylation, breaks=c(0, 0.0000001, 0.25, 0.5, 0.75, 1), include.lowest=T)))
DKO_0.35uM_1i_Epi <- data.frame(table(cut(DKO_0.35uM_1i_Epi_CRs$mean_methylation, breaks=c(0, 0.0000001, 0.25, 0.5, 0.75, 1), include.lowest=T)))

Dnmt1_KO_E3.5$frac <- 100*Dnmt1_KO_E3.5$Freq/sum(Dnmt1_KO_E3.5$Freq)
Dnmt1_KO_E3.5$type <- 'E3.5'
Dnmt1_KO_E3.5$group <- 'Dnmt1_KO'
Dnmt1_KO_E6.5_Epi$frac <- 100*Dnmt1_KO_E6.5_Epi$Freq/sum(Dnmt1_KO_E6.5_Epi$Freq)
Dnmt1_KO_E6.5_Epi$type <- 'E6.5_Epi'
Dnmt1_KO_E6.5_Epi$group <- 'Dnmt1_KO'
DKO_0.35uM_1i_E35$frac <- 100*DKO_0.35uM_1i_E35$Freq/sum(DKO_0.35uM_1i_E35$Freq)
DKO_0.35uM_1i_E35$type <- 'E3.5'
DKO_0.35uM_1i_E35$group <- 'DKO_0.35uM_1i'
DKO_0.35uM_1i_Epi$frac <- 100*DKO_0.35uM_1i_Epi$Freq/sum(DKO_0.35uM_1i_Epi$Freq)
DKO_0.35uM_1i_Epi$type <- 'E6.5_Epi'
DKO_0.35uM_1i_Epi$group <- 'DKO_0.35uM_1i'

data <- rbind(Dnmt1_KO_E3.5, Dnmt1_KO_E6.5_Epi, DKO_0.35uM_1i_E35, DKO_0.35uM_1i_Epi)
data$type <- factor(data$type, levels=c('E3.5', 'E6.5_Epi'))
data$Var1 <- factor(data$Var1, levels=rev(levels(data$Var1)))

pdf('read_methylation_level_in_vivo_CRs.pdf')
ggplot(data, aes(x = type, stratum = Var1, alluvium = Var1, y = frac, fill = Var1, label = Var1)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow(alpha = .5) +
  geom_stratum(color = 'lightgrey') +
  theme_classic() +
  scale_fill_manual(values=rev(c('dodgerblue3','chartreuse3','yellow','orange','firebrick'))) + 
  facet_grid(~group)
dev.off()
```
