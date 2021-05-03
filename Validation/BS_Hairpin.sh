```R
require(ggplot2)
require(reshape2)

data <- read.table('BS_Hairpin.tsv', header=T)
data <- dcast(melt(data), Amplicon+variable~Sample)

data$DKO_P1 <- log2(data[,3]/data[,5])
data$DKO_P5 <- log2(data[,4]/data[,5])
data$WT <- log2(data[,6]/data[,5])

df <- melt(data[,c(1,2,7,8,9)])
colnames(df) <- c('Amplicon','Methylation','Sample','log2FC_percent')

pdf('figures/hairpin_quant.pdf', width=4)
ggplot(subset(df, Methylation != 'fully_methylated' & Amplicon %in% c('iapez','iapltr','l1mdt','msat','muervl')), aes(x=Methylation, y=log2FC_percent, fill=Sample)) + geom_bar(stat='identity', color='black', position='dodge') + theme_classic() + facet_grid(Amplicon~.)
dev.off()
```
