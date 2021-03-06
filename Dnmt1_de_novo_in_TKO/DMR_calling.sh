# Characterization methylation differences genome-wide
for i in TKOlike.bw DKOzero_P15.bw; do
	n=`basename $i | sed 's/.bw//'`
	bigWigAverageOverBed ${i} w1000.bed ${n}_1kb.tsv
done

```R
require(ggplot2)
require(reshape2)

# load data
TKO <- read.table('TKOlike_1kb.tsv', header=F, col.names=c('name','size','covered','sum','mean0','TKO'))
DKO <- read.table('/DKOzero_P15_1kb.tsv', header=F, col.names=c('name','size','covered','sum','mean0','DKO'))

df <- merge(subset(TKO, covered>=5)[,c(1,6)], subset(DKO, covered>=5)[,c(1,6)], by='name')
df$d <- df$DKO-df$TKO

pdf('figures/TKO_DKO_1kb_methylation_density.pdf', width=14)
ggplot(melt(df[,c(1:3)]), aes(x=value, color=variable)) + geom_line(stat='density', adjust=0.05) + theme_classic() + facet_wrap(~variable, scales='free_y') + xlim(c(0,0.5)) + geom_vline(xintercept=c(0.05))
dev.off()

pdf('figures/TKO_DKO_1kb_methylation_difference_ecf.pdf')
ggplot(df, aes(d)) + stat_ecdf(size=2) + theme_classic() + geom_hline(yintercept=0.95, color='grey') + geom_vline(xintercept=0.1, color='grey') + xlim(c(0,0.5))
dev.off()
```


# DMR calling

## sliding window w 1000 s 250
### average methylation
bedtools intersect -wa -wb -a w1000_s250.bed -b DKO_TKO_diff.bed | cut -f1,2,3,4,8 | bedtools groupby -g 1,2,3,4 -c 5 -o median,max,count,collapse >diff_w1000_s250_median.txt

### differential regions
```R
require(ggplot2)

diff <- read.table('diff_w1000_s250_median.txt', header=F, col.names=c('chr','start','end','w','me','max','count','rates'))
diff <- subset(diff, count>=10)

pdf('figures/w1000_s250_median.pdf')
ggplot(diff, aes(count)) + geom_histogram(bins=50) + theme_classic() + scale_x_continuous(limits=c(0,50)) + ggtitle('DKO-TKO') + xlab('#CpGs')
ggplot(diff, aes(me)) + geom_histogram(bins=100) + theme_classic() + scale_x_continuous(limits=c(0,1)) + ggtitle('DKO-TKO') + xlab('Avg methylation')
dev.off()

write.table(file='diff.0.15_w1000_s250_median.bed', subset(diff, abs(me) > 0.15)[,c('chr','start','end','me')], sep='\t', quote=F, col.names=F, row.names=F)
dim(subset(diff, abs(me) > 0.15))
```

## sliding window w 2000 s 500
### average methylation
bedtools intersect -wa -wb -a w2000_s500.bed -b DKO_TKO_diff.bed | cut -f1,2,3,4,8 | bedtools groupby -g 1,2,3,4 -c 5 -o median,max,count,collapse >diff_w2000_s500_median.txt

### differential regions
```R
require(ggplot2)

diff <- read.table('diff_w2000_s500_median.txt', header=F, col.names=c('chr','start','end','w','me','max','count','rates'))
diff <- subset(diff, count>=10)

pdf('figures/w2000_s500_median.pdf')
ggplot(diff, aes(count)) + geom_histogram(bins=50) + theme_classic() + scale_x_continuous(limits=c(0,50)) + ggtitle('DKO-TKO') + xlab('#CpGs')
ggplot(diff, aes(me)) + geom_histogram(bins=100) + theme_classic() + scale_x_continuous(limits=c(0,1)) + ggtitle('DKO-TKO') + xlab('Avg methylation')
dev.off()

write.table(file='diff.0.15_w2000_s500_median.bed', subset(diff, abs(me) > 0.15)[,c('chr','start','end','me')], sep='\t', quote=F, col.names=F, row.names=F)
dim(subset(diff, abs(me) > 0.15))
```

## merge peaks and filter for DMRs with an average methylation in P15 > 0.15
less diff.0.15_w1000_s250_median.bed diff.0.15_w2000_s500_median.bed | bedtools sort | bedtools merge -d 5000 -c 4 -o mean | bedtools intersect -v -a stdin -b Nanopore_deletions.bed >tmp
bedtools intersect -wa -wb -a tmp -b DKOzero_P15.bed | bedtools groupby -g 1,2,3,4 -c 8 -o mean | perl -ane 'if($F[4]>0.15){print $_}' | cut -f1-4 >DKO_TKO_DMRs.bed

less DKO_TKO_DMRs.bed | perl -ane 'BEGIN{$i=1} print "$F[0]\t$F[1]\t$F[2]\tDMR_$i;$F[3]\n"; $i++' >DKO_TKO_DMRs_numerated_v1.bed


# CRs
## randomly shuffle DMRs for control set
bedtools merge -d 2500 -i DKO_TKO_diff.bed -c 4 -o count | perl -ane 'if($F[3]>=10){print $_}' >keep
for r in {1..1000}; do
    echo "Iteration $r" >&2
    bedtools shuffle -i DKO_TKO_DMRs_numerated_v1.bed -g mm9.chrom.sizes -incl keep | awk -v r=$r '{print $0";iteration_"r}'
done >DKO_TKO_CRs_numerated_v1.bed


# DMR and CR length and differences
bedtools intersect -wao -a DKO_TKO_DMRs_numerated_v1.bed -b DKO_TKO_diff.bed | bedtools groupby -g 1,2,3,4 -c 8,9 -o mean,count | perl -ane 'BEGIN{print "#chr\tstart\tend\tID\tdiff\tCpG\tlength\n"}; $l=$F[2]-$F[1]; print "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$l\n"' >DKO_TKO_DMRs_length.txt
bedtools intersect -wao -a DKO_TKO_CRs_numerated_v1.bed -b DKO_TKO_diff.bed | bedtools groupby -g 1,2,3,4 -c 8,9 -o mean,count | perl -ane 'BEGIN{print "#chr\tstart\tend\tID\tdiff\tCpG\tlength\n"}; $l=$F[2]-$F[1]; print "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$l\n"' >DKO_TKO_CRs_length.txt
bedtools intersect -wao -a cr_matching_cpg_num/DKO_TKO_CRs_CpGdensity_v1.bed -b DKO_TKO_diff.bed | bedtools groupby -g 1,2,3,4 -c 8,9 -o mean,count | perl -ane 'BEGIN{print "#chr\tstart\tend\tID\tdiff\tCpG\tlength\n"}; $l=$F[2]-$F[1]; print "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$l\n"' >DKO_TKO_CRs2_length.txt

bedtools intersect -wao -a w2000_s500.bed -b DKO_TKO_diff.bed | bedtools groupby -g 1,2,3,4 -c 8,9 -o mean,count | perl -ane 'BEGIN{print "#chr\tstart\tend\tID\tdiff\tCpG\tlength\n"}; $l=$F[2]-$F[1]; print "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$l\n"' >DKO_TKO_w_length.txt


```R
require(ggplot2)
require(ggExtra)

DMRs <- read.table('DKO_TKO_DMRs_length.txt', header=T, comment.char='')
CRs <- read.table('DKO_TKO_CRs_length.txt', header=T, comment.char='')
CRs2 <- read.table('DKO_TKO_CRs2_length.txt', header=T, comment.char='')

p <- ggplot(DMRs, aes(x=length, y=CpG, color=diff)) + geom_point() + theme_classic() + xlim(c(0,22000)) + ylim(c(0,450)) + scale_color_gradientn(colours = c('darkblue','gray90','firebrick3'), limits=c(-0.35,0.35), space='Lab') + xlab('DMR length [nt]') + ylab('DMR length [CpGs]') + geom_abline(intercept=0, slope=c(0.005, 0.008))

pdf('figures/DKO_TKO_DMRs_length.pdf')
ggMarginal(p, size=4, xparams = list(size = 1), yparams = list(size = 1))
plot.new()
ggMarginal(p + theme(legend.position='none'), size=4, xparams = list(size = 1), yparams = list(size = 1))
dev.off()


set.seed(12)
x <- table(DMRs$length)
subCRs <- data.frame(X.chr=factor(), start=integer(), end=integer(), ID=integer(), diff=numeric(), CpG=integer(), length=integer())
for (i in 1:length(x)){
	subCRs <- rbind(subCRs, CRs[sample(which(CRs$length == as.numeric(names(x[i]))), x[i]), ])
}

p <- ggplot(subCRs, aes(x=length, y=CpG, color=diff)) + geom_point() + theme_classic() + xlim(c(0,22000)) + ylim(c(0,450)) + scale_color_gradientn(colours = c('darkblue','gray90','firebrick3'), limits=c(-0.35,0.35), space='Lab') + xlab('CR length [nt]') + ylab('CR length [CpGs]') + geom_abline(intercept=0, slope=c(0.005, 0.008))

pdf('figures/DKO_TKO_CRs_length.pdf')
ggMarginal(p, size=4, xparams = list(size = 1), yparams = list(size = 1))
plot.new()
ggMarginal(p + theme(legend.position='none'), size=4, xparams = list(size = 1), yparams = list(size = 1))
dev.off()


set.seed(12)
x <- table(DMRs$length)
subCRs <- data.frame(X.chr=factor(), start=integer(), end=integer(), ID=integer(), diff=numeric(), CpG=integer(), length=integer())
for (i in 1:length(x)){
	subCRs <- rbind(subCRs, CRs2[sample(which(CRs2$length == as.numeric(names(x[i]))), x[i]), ])
}

p <- ggplot(subCRs, aes(x=length, y=CpG, color=diff)) + geom_point() + theme_classic() + xlim(c(0,22000)) + ylim(c(0,450)) + scale_color_gradientn(colours = c('darkblue','gray90','firebrick3'), limits=c(-0.35,0.35), space='Lab') + xlab('CR length [nt]') + ylab('CR length [CpGs]') + geom_abline(intercept=0, slope=c(0.005, 0.008))

pdf('figures/DKO_TKO_CRs2_length.pdf')
ggMarginal(p, size=4, xparams = list(size = 1), yparams = list(size = 1))
plot.new()
ggMarginal(p + theme(legend.position='none'), size=4, xparams = list(size = 1), yparams = list(size = 1))
dev.off()
```
