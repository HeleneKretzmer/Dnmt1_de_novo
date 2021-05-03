# libraries
require(ggplot2)
require(plyr)
require(reshape2)


# load data
Nuclei <- read.csv("Nuclei.csv", sep=";", comment.char = '#', dec=',')

# add meta info
Nuclei$sample <- gsub('^1$','DKOP25',gsub('.*_','',gsub('-Orthogonal.*','',Nuclei$ImageDocumentName..Image.Name)))
Nuclei$sample <- factor(Nuclei$sample, levels=c('WT','TKOL','DKOP5','DKOP15','DKOP25'))
#Nuclei$projection <- gsub('.czi','',gsub('.*-Orthogonal ','',Nuclei$ImageDocumentName..Image.Name))

# replace NA for measurements by 0
Nuclei[is.na(Nuclei)] <- 0

# plots
ggplot(Nuclei, aes(x=sample, y=log10(Cy5_Spots.RegionsArea...Cy5_Spots.Area..R))) + geom_boxplot() + theme_classic()
ggplot(Nuclei, aes(x=sample, y=Cy5_Spots.RegionsCount...Cy5_Spots.Count..I)) + geom_boxplot() + theme_classic()


df <- ddply(Nuclei, .(sample, ImageSceneName..Image.Scene.Name), summarize, m=mean((IntensityMean_Cy5..Intensity.Mean.Value.of.channel..Cy5...R)/ImageSceneName..Image.Scene.Name))
ggplot(df, aes(x=sample, y=log10(m))) + geom_boxplot() + theme_classic()

df <- ddply(Nuclei, .(sample, ImageSceneName..Image.Scene.Name), summarize, m=mean(Area..Area..R))
ggplot(df, aes(x=sample, y=m)) + geom_boxplot() + theme_classic()


# randomly sample 100 cells/measurement
set.seed(12)
n <- 100

IAPpos <- data.frame(WT=numeric(), TKOL=numeric(), DKOP5=numeric(), DKOP15=numeric(), DKOP25=numeric())
for (i in 1:1000){
	WTcells <- sample(which(Nuclei$sample=='WT'), n)
	TKOLcells <- sample(which(Nuclei$sample=='TKOL'), n)
	P5cells <- sample(which(Nuclei$sample=='DKOP5'), n)
	P15cells <- sample(which(Nuclei$sample=='DKOP15'), n)
	P25cells <- sample(which(Nuclei$sample=='DKOP25'), n)

	# IAP pos cells
	WTcells_IAP <- sum(Nuclei[WTcells,]$Cy5_Spots.RegionsCount...Cy5_Spots.Count..I > 0)
	TKOLcells_IAP <- sum(Nuclei[TKOLcells,]$Cy5_Spots.RegionsCount...Cy5_Spots.Count..I > 0)
	P5cells_IAP <- sum(Nuclei[P5cells,]$Cy5_Spots.RegionsCount...Cy5_Spots.Count..I > 0)
	P15cells_IAP <- sum(Nuclei[P15cells,]$Cy5_Spots.RegionsCount...Cy5_Spots.Count..I > 0)
	P25cells_IAP <- sum(Nuclei[P25cells,]$Cy5_Spots.RegionsCount...Cy5_Spots.Count..I > 0)
	IAPpos <- rbind(IAPpos, data.frame(WT=WTcells_IAP, TKOL=TKOLcells_IAP, DKOP5=P5cells_IAP, DKOP15=P15cells_IAP, DKOP25=P25cells_IAP))
}


# plotting
ggplot(melt(IAPpos), aes(x=variable, y=value)) + geom_boxplot() + theme_classic()
