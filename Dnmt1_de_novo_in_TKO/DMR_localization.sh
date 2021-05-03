# Genome wide distribution (http://liulab.dfci.harvard.edu/CEAS/usermanual.html)
source virtual/CEAS/bin/activate
cat DKO_TKO_DMRs_numerated_v1.bed | sed 's/;/\t/' >DKO_TKO_DMRs_genome.annotation
ceas -g mm9.refGene -b DKO_TKO_DMRs_genome.annotation
deactivate


# Overlap with genomic elements
bedtools merge -i CGI.bed | bedtools intersect -wao -a DKO_TKO_DMRs.bed -b stdin | bedtools groupby -g 1,2,3,4 -c 8 -o sum | perl -ane 'BEGIN{print "DMR\tregion\tlength\toverlap\n"}; $l=$F[2]-$F[1];  print "DMR\tCGI\t$l\t$F[-1]\n"' >DKO_TKO_DMRs_CGI.txt
bedtools sort -i gencode.protein_coding.TSS.bed | bedtools merge | bedtools intersect -wao -a DKO_TKO_DMRs.bed -b stdin | bedtools groupby -g 1,2,3,4 -c 8 -o sum | perl -ane 'BEGIN{print "DMR\tregion\tlength\toverlap\n"}; $l=$F[2]-$F[1];  print "DMR\tTSS\t$l\t$F[-1]\n"' >DKO_TKO_DMRs_TSS.txt
bedtools merge -i LINE.bed | bedtools intersect -wao -a DKO_TKO_DMRs.bed -b stdin | bedtools groupby -g 1,2,3,4 -c 8 -o sum | perl -ane 'BEGIN{print "DMR\tregion\tlength\toverlap\n"}; $l=$F[2]-$F[1];  print "DMR\tLINE\t$l\t$F[-1]\n"' >DKO_TKO_DMRs_LINE.txt
bedtools merge -i SINE.bed | bedtools intersect -wao -a DKO_TKO_DMRs.bed -b stdin | bedtools groupby -g 1,2,3,4 -c 8 -o sum | perl -ane 'BEGIN{print "DMR\tregion\tlength\toverlap\n"}; $l=$F[2]-$F[1];  print "DMR\tSINE\t$l\t$F[-1]\n"' >DKO_TKO_DMRs_SINE.txt
bedtools merge -i LTR.bed | bedtools intersect -wao -a DKO_TKO_DMRs.bed -b stdin | bedtools groupby -g 1,2,3,4 -c 8 -o sum | perl -ane 'BEGIN{print "DMR\tregion\tlength\toverlap\n"}; $l=$F[2]-$F[1];  print "DMR\tLTR\t$l\t$F[-1]\n"' >DKO_TKO_DMRs_LTR.txt

bedtools merge -i CGI.bed | bedtools intersect -wao -a DKO_TKO_CRs.bed -b stdin | bedtools groupby -g 1,2,3,4 -c 8 -o sum | perl -ane 'BEGIN{print "DMR\tregion\tlength\toverlap\n"}; $l=$F[2]-$F[1];  print "CR\tCGI\t$l\t$F[-1]\n"' >DKO_TKO_CRs_CGI.txt
bedtools sort -i gencode.protein_coding.TSS.bed | bedtools merge | bedtools intersect -wao -a DKO_TKO_CRs.bed -b stdin | bedtools groupby -g 1,2,3,4 -c 8 -o sum | perl -ane 'BEGIN{print "DMR\tregion\tlength\toverlap\n"}; $l=$F[2]-$F[1];  print "CR\tTSS\t$l\t$F[-1]\n"' >DKO_TKO_CRs_TSS.txt
bedtools merge -i LINE.bed | bedtools intersect -wao -a DKO_TKO_CRs.bed -b stdin | bedtools groupby -g 1,2,3,4 -c 8 -o sum | perl -ane 'BEGIN{print "DMR\tregion\tlength\toverlap\n"}; $l=$F[2]-$F[1];  print "CR\tLINE\t$l\t$F[-1]\n"' >DKO_TKO_CRs_LINE.txt
bedtools merge -i SINE.bed | bedtools intersect -wao -a DKO_TKO_CRs.bed -b stdin  | bedtools groupby -g 1,2,3,4 -c 8 -o sum| perl -ane 'BEGIN{print "DMR\tregion\tlength\toverlap\n"}; $l=$F[2]-$F[1];  print "CR\tSINE\t$l\t$F[-1]\n"' >DKO_TKO_CRs_SINE.txt
bedtools merge -i LTR.bed | bedtools intersect -wao -a DKO_TKO_CRs.bed -b stdin | bedtools groupby -g 1,2,3,4 -c 8 -o sum | perl -ane 'BEGIN{print "DMR\tregion\tlength\toverlap\n"}; $l=$F[2]-$F[1];  print "CR\tLTR\t$l\t$F[-1]\n"' >DKO_TKO_CRs_LTR.txt

bedtools merge -i CGI.bed | bedtools intersect -wao -a DKO_TKO_CRs_CpGdensity_v1.bed -b stdin | bedtools groupby -g 1,2,3,4 -c 8 -o sum | perl -ane 'BEGIN{print "DMR\tregion\tlength\toverlap\n"}; $l=$F[2]-$F[1];  print "CR2\tCGI\t$l\t$F[-1]\n"' >DKO_TKO_CRs2_CGI.txt
bedtools sort -i gencode.protein_coding.TSS.bed | bedtools merge | bedtools intersect -wao -a DKO_TKO_CRs_CpGdensity_v1.bed -b stdin | bedtools groupby -g 1,2,3,4 -c 8 -o sum | perl -ane 'BEGIN{print "DMR\tregion\tlength\toverlap\n"}; $l=$F[2]-$F[1];  print "CR2\tTSS\t$l\t$F[-1]\n"' >DKO_TKO_CRs2_TSS.txt
bedtools merge -i LINE.bed | bedtools intersect -wao -a DKO_TKO_CRs_CpGdensity_v1.bed -b stdin | bedtools groupby -g 1,2,3,4 -c 8 -o sum | perl -ane 'BEGIN{print "DMR\tregion\tlength\toverlap\n"}; $l=$F[2]-$F[1];  print "CR2\tLINE\t$l\t$F[-1]\n"' >DKO_TKO_CRs2_LINE.txt
bedtools merge -i SINE.bed | bedtools intersect -wao -a DKO_TKO_CRs_CpGdensity_v1.bed -b stdin | bedtools groupby -g 1,2,3,4 -c 8 -o sum| perl -ane 'BEGIN{print "DMR\tregion\tlength\toverlap\n"}; $l=$F[2]-$F[1];  print "CR2\tSINE\t$l\t$F[-1]\n"' >DKO_TKO_CRs2_SINE.txt
bedtools merge -i LTR.bed | bedtools intersect -wao -a DKO_TKO_CRs_CpGdensity_v1.bed -b stdin | bedtools groupby -g 1,2,3,4 -c 8 -o sum | perl -ane 'BEGIN{print "DMR\tregion\tlength\toverlap\n"}; $l=$F[2]-$F[1];  print "CR2\tLTR\t$l\t$F[-1]\n"' >DKO_TKO_CRs2_LTR.txt

```R
require(plyr)
require(ggplot2)

DMR_CGI <- read.table('DKO_TKO_DMRs_CGI.txt', header=T)
CR_CGI <- read.table('DKO_TKO_CRs_CGI.txt', header=T)
CR2_CGI <- read.table('DKO_TKO_CRs2_CGI.txt', header=T)
CGI_pval <- wilcox.test(DMR_CGI[,4], CR_CGI[,4], alternative='greater')$p.val

DMR_TSS <- read.table('DKO_TKO_DMRs_TSS.txt', header=T)
CR_TSS <- read.table('DKO_TKO_CRs_TSS.txt', header=T)
CR2_TSS <- read.table('DKO_TKO_CRs2_TSS.txt', header=T)
TSS_pval <- wilcox.test(DMR_TSS[,4], CR_TSS[,4], alternative='greater')$p.val

DMR_LINE <- read.table('DKO_TKO_DMRs_LINE.txt', header=T)
CR_LINE <- read.table('DKO_TKO_CRs_LINE.txt', header=T)
CR2_LINE <- read.table('DKO_TKO_CRs2_LINE.txt', header=T)
LINE_pval <- wilcox.test(DMR_LINE[,4], CR_LINE[,4], alternative='greater')$p.val

DMR_SINE <- read.table('DKO_TKO_DMRs_SINE.txt', header=T)
CR_SINE <- read.table('DKO_TKO_CRs_SINE.txt', header=T)
CR2_SINE <- read.table('DKO_TKO_CRs2_SINE.txt', header=T)
SINE_pval <- wilcox.test(DMR_SINE[,4], CR_SINE[,4], alternative='greater')$p.val

DMR_LTR <- read.table('DKO_TKO_DMRs_LTR.txt', header=T)
CR_LTR <- read.table('DKO_TKO_CRs_LTR.txt', header=T)
CR2_LTR <- read.table('DKO_TKO_CRs2_LTR.txt', header=T)
LTR_pval <- wilcox.test(DMR_LTR[,4], CR_LTR[,4], alternative='greater')$p.val


df <- rbind(DMR_TSS,CR_TSS, DMR_CGI,CR_CGI, DMR_LINE,CR_LINE, DMR_SINE,CR_SINE, DMR_LTR,CR_LTR)
df <- ddply(df, .(DMR, region), summarize, frac=100*sum(overlap)/sum(length))

pdf('figures/DKO_TKO_genomic_regions_enrichment.pdf')
ggplot(df, aes(x=region, fill=DMR, y=frac)) + geom_bar(stat='identity', position='dodge', color='black') + theme_classic() + scale_fill_manual(values=c('gold1','grey')) + ylab('Fraction of nt overlap')
dev.off()

stats <- data.frame(Region=c('CGI','TSS','LINE','SINE','LTR'), pval=c(CGI_pval,TSS_pval,LINE_pval,SINE_pval,LTR_pval))


df <- rbind(DMR_TSS,CR_TSS,CR2_TSS, DMR_CGI,CR_CGI,CR2_CGI, DMR_LINE,CR_LINE,CR2_LINE, DMR_SINE,CR_SINE,CR2_SINE, DMR_LTR,CR_LTR,CR2_LTR)
df <- ddply(df, .(DMR, region), summarize, frac=100*sum(overlap)/sum(length))

pdf('figures/DKO_TKO_genomic_regions_enrichment_v2.pdf')
ggplot(df, aes(x=region, fill=DMR, y=frac)) + geom_bar(stat='identity', position='dodge', color='black') + theme_classic() + scale_fill_manual(values=c('gold1','grey','darkgrey')) + ylab('Fraction of nt overlap')
dev.off()
```

# Overlap with LTR classes
echo -ne "region\tDMR\tlength\tchr\tstart\tend\tID\tClass\tFamily\toverlap\n" >DKO_TKO_DMRs_LTR_Class_nt.txt
echo -ne "region\tDMR\tlength\tchr\tstart\tend\tID\tClass\tFamily\toverlap\n" >DKO_TKO_CRs_LTR_Class_nt.txt
for i in `cat LTR.bed | cut -f5 | sort | uniq`; do
    echo $i
    cat LTR.bed | awk -v r=$i '{ if ($5 == r) print $0 }' | bedtools intersect -wao -a DKO_TKO_DMRs.bed -b stdin | cut -f1-4,9,10,13 | bedtools groupby -g 1,2,3,4,5,6 -c 7 -o sum | perl -ane '$l=$F[2]-$F[1]; print "DMR\t$l\t$_"' | awk -v r=$i '{print r"\t"$0}' >>DKO_TKO_DMRs_LTR_Class_nt.txt
    cat LTR.bed | awk -v r=$i '{ if ($5 == r) print $0 }' | bedtools intersect -wao -a DKO_TKO_CRs.bed -b stdin | cut -f1-4,9,10,13 | bedtools groupby -g 1,2,3,4,5,6 -c 7 -o sum | perl -ane '$l=$F[2]-$F[1]; print "CR\t$l\t$_"' | awk -v r=$i '{print r"\t"$0}' >>DKO_TKO_CRs_LTR_Class_nt.txt
done

```R
require(plyr)
require(ggplot2)

DMR_LTR <- read.table('DKO_TKO_DMRs_LTR_Class_nt.txt', header=T)
CR_LTR <- read.table('DKO_TKO_CRs_LTR_Class_nt.txt', header=T)

DMR <- ddply(DMR_LTR, .(DMR, Class), summarize, nt=sum(overlap))
CR <- ddply(CR_LTR, .(DMR, Class), summarize, nt=sum(overlap))

l_DMR <- sum(unique(DMR_LTR[,c(4:7)])$end - unique(DMR_LTR[,c(4:7)])$start)
l_CR <- sum(unique(CR_LTR[,c(4:7)])$end - unique(CR_LTR[,c(4:7)])$start)
DMR$frac <- 100*DMR$nt/l_DMR
CR$frac <- 100*CR$nt/l_CR

df <- rbind(DMR, CR)

pdf('figures/DKO_TKO_LTR_Class_enrichment.pdf', width=14)
ggplot(df, aes(x=Class, fill=DMR, y=frac)) + geom_bar(stat='identity', position='dodge', color='black') + theme_classic() + scale_fill_manual(values=c('gold1','grey')) + ylab('Fraction of nt overlap') + ylim(c(0,80))
dev.off()

for (r in levels(df$Class)[-1]){
    cat(r, '\n')
    pval <- wilcox.test(subset(DMR_LTR, region == r)$overlap, subset(CR_LTR, region == r)$overlap, alternative='greater')$p.val
    assign(paste0(r, '_pval'), pval)
}

stats <- data.frame(Region=c('ERV1', 'ERVK', 'ERVL', 'Gypsy', 'MaLR', 'ERV', 'LTR'), pval=c(ERV1_pval, ERVK_pval, ERVL_pval, Gypsy_pval, MaLR_pval, ERV_pval, LTR_pval))
stats$nt <- df$nt[match(stats$Region, df$Class)]
stats$frac <- df$frac[match(stats$Region, df$Class)]
```


## IAPs vs all other ERVK
echo -ne "region\tDMR\tClass\tlength\tchr\tstart\tend\tID\tClass\tFamily\toverlap\n" >DKO_TKO_DMRs_ERVK_nt.txt
echo -ne "region\tDMR\tClass\tlength\tchr\tstart\tend\tID\tClass\tFamily\toverlap\n" >DKO_TKO_CRs_ERVK_nt.txt
for i in `cat ERVK.bed | cut -f4 | sort | uniq`; do
    echo $i
    cat ERVK.bed | awk -v r=$i '{ if ($4 == r) print $0 }' | bedtools intersect -wao -a DKO_TKO_DMRs.bed -b stdin | cut -f1-4,8,9,13 | bedtools groupby -g 1,2,3,4,5,6 -c 7 -o sum | perl -ane '$l=$F[2]-$F[1]; print "DMR\tLTR\t$l\t$_"' | awk -v r=$i '{print r"\t"$0}' >>DKO_TKO_DMRs_ERVK_nt.txt
    cat ERVK.bed | awk -v r=$i '{ if ($4 == r) print $0 }' | bedtools intersect -wao -a DKO_TKO_CRs.bed -b stdin | cut -f1-4,8,9,13 | bedtools groupby -g 1,2,3,4,5,6 -c 7 -o sum | perl -ane '$l=$F[2]-$F[1]; print "CR\tLTR\t$l\t$_"' | awk -v r=$i '{print r"\t"$0}' >>DKO_TKO_CRs_ERVK_nt.txt
done

```R
require(plyr)
require(ggplot2)

DMR_LTR <- read.table('DKO_TKO_DMRs_ERVK_nt.txt', header=T, nrow=227251)
CR_LTR <- read.table('DKO_TKO_CRs_ERVK_nt.txt', header=T, nrow=227250001)

DMR <- ddply(DMR_LTR, .(DMR, Class.1), summarize, nt=sum(overlap))
CR <- ddply(CR_LTR, .(DMR, Class.1), summarize, nt=sum(overlap))

l_DMR <- sum(unique(DMR_LTR[,c(4:7)])$end - unique(DMR_LTR[,c(4:7)])$start)
l_CR <- 1000*l_DMR
DMR$frac <- 100*DMR$nt/l_DMR
CR$frac <- 100*CR$nt/l_CR

df <- rbind(DMR, CR)

pdf('figures/DKO_TKO_ERVK_enrichment.pdf', width=14)
ggplot(subset(df, Class.1 %in% subset(df, frac>0.1)$Class.1), aes(x=Class.1, fill=DMR, y=frac)) + geom_bar(stat='identity', position='dodge', color='black') + theme_classic() + scale_fill_manual(values=c('gold1','grey')) + ylab('Fraction of nt overlap') + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + xlab('ERVK elements overlapping at least 10% of DMRs/CRs')
dev.off()

for (r in unique(subset(df, Class.1 %in% subset(df, frac>0.1)$Class.1)$Class.1)){
    cat(r, '\n')
    pval <- wilcox.test(subset(DMR_LTR, region == r)$overlap, subset(CR_LTR, region == r)$overlap, alternative='greater')$p.val
    assign(paste0(r, '_pval'), pval)
    cat(pval, '\n')
}

stats <- data.frame(Region=c('IAP-d-int', 'IAPEY2_LTR', 'IAPEY3-int', 'IAPEy-int', 'IAPEz-int', 'IAPLTR1_Mm', 'IAPLTR1a_Mm', 'IAPLTR2_Mm', 'IAPLTR2a', 'IAPLTR2b', 'IAPLTR3-int', 'IAPLTR4_I', 'MMERVK10C-int', 'MMETn-int', 'RLTR10', 'RLTR10-int', 'RLTR10C', 'RLTR10D', 'RLTR27', 'RLTR44-int', 'RLTR45-int', 'RLTRETN_Mm', 'RMER16-int', 'RMER17B', 'RMER19C', 'RMER6C'), pval=c(IAP-d-int_pval, IAPEY2_LTR_pval, IAPEY3-int_pval, IAPEy-int_pval, IAPEz-int_pval, IAPLTR1_Mm_pval, IAPLTR1a_Mm_pval, IAPLTR2_Mm_pval, IAPLTR2a_pval, IAPLTR2b_pval, IAPLTR3-int_pval, IAPLTR4_I_pval, MMERVK10C-int_pval, MMETn-int_pval, RLTR10_pval, RLTR10-int_pval, RLTR10C_pval, RLTR10D_pval, RLTR27_pval, RLTR44-int_pval, RLTR45-int_pval, RLTRETN_Mm_pval, RMER16-int_pval, RMER17B_pval, RMER19C_pval, RMER6C_pval))
```
