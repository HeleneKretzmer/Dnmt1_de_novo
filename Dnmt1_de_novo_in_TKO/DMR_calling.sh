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


