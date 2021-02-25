#### Code by AEB ####
#start at Line 136 for QTL analysis, data for this provided

##### genetic map generation for RILs pop ####
library(qtl)

#The genotypes of this mapping population can only be released upon direct request
#please email cris.kuhlemeier@ips.unibe.ch for data requests

mapthis <- read.cross("csv", "~/RILs/maps/", "Peax.mathieu.perf.out.Andrea.physmap.csv", estimate.map = F, crosstype = "riself")

plotMissing(mapthis)

#remove individuals with very few genotypes
mapthis <- subset(mapthis, ind=(ntyped(mapthis)>75)) #remove individuals with a lot of missing
summary(mapthis)
#tot 192 ind left

#To omit the markers with lots of missing data, we first need to identify the names of themarkers
nt.bymar <- ntyped(mapthis, "mar")
todrop <- names(nt.bymar[nt.bymar < 100])
mapthis <- drop.markers(mapthis, todrop)
summary(mapthis)
#none dropped

#identify duplicate individuals
cg <- comparegeno(mapthis)
hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes")
rug(cg[lower.tri(cg)])

wh <- which(cg > 0.9, arr=TRUE) #identify pairs that have over 90% matching genos
wh <- wh[wh[,1] < wh[,2],]
wh
#eh, that's a lot, let's keep them in. This isn't a perfect mapping pop.

#duplicate markers
print(dup <- findDupMarkers(mapthis, exact.only=FALSE))
#lots, keep, not a perfect mapping pop.

#markers with seg dist
gt <- geno.table(mapthis)
gt[gt$P.value < 0.05/totmar(mapthis),]
#there's a lot, but we knew this. drop the worst ones
todrop <- rownames(gt[gt$P.value < 1e-10,])
mapthis <- drop.markers(mapthis, todrop)
summary(mapthis)

#individuals' genotype frequencies, modified for RIself (no hets)
g <- pull.geno(mapthis)
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:3)))
gfreq <- t(t(gfreq) / colSums(gfreq))
par(mfrow=c(1,3), las=1)
for(i in 1:3)
  plot(gfreq[i,], ylab="Genotype frequency", main=c("AA", "BB")[i], ylim=c(0,1))
par(mfrow=c(1,1))
#Study pairwise marker linkages; look for switched alleles
mapthis <- est.rf(mapthis)
checkAlleles(mapthis, threshold=5)
#no apparent problems
#no switched alleles
rf <- pull.rf(mapthis)
lod <- pull.rf(mapthis, what="lod")
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")

#let's try max.rf-0.35, min.lod=6
lg <- formLinkageGroups(mapthis, max.rf=0.35, min.lod=6)
table(lg[,2])
#  1   2   3   4   5   6   7   8   9  10 
#279 279 250 214 184 151  27  19   5   1 
#pretty good
#tried other combinations that weren't so good, keeping this

lg <- formLinkageGroups(mapthis, max.rf=0.35, min.lod=5)
table(lg[,2])
# 1   2   3   4   5   6   7   8   9 
# 279 279 250 214 184 152  27  19   5 
#tried other combinations that weren't so good, keeping this

mapthis <- formLinkageGroups(mapthis, max.rf = 0.35, min.lod = 5, reorgMarkers = T)
plotRF(mapthis, alternate.chrid=TRUE)


mapthis <- formLinkageGroups(mapthis, max.rf = 0.35, min.lod = 6, reorgMarkers = T)
plotRF(mapthis, alternate.chrid=TRUE)
#good enough for now

#what chromosomes are which LGs
mapthis$geno$`1`
chrnames(mapthis)
names(mapthis$geno)
#[1] "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10"
names(mapthis$geno) <- c("2", "3", "4", "1", "7", "5", "6", "8", "9", "10")
    #renamed chromosomes based on marker locations, as noted in their names

#LG8 has chr 1 markers
#drop LG9 and LG10
plotRF(mapthis, alternate.chrid=TRUE)


library(ASMap)
#using ASMap to order the markers in LGs - it is much faster

orderedmap <- mstmap(mapthis, id="Line", bychr = TRUE, dist.fun = "kosambi", anchor = FALSE,
                     trace = TRUE)
nmar(orderedmap)
chrlen(orderedmap)
heatMap(orderedmap, lmax = 70)



#To ensure the extra list elements "co.located", "seg.distortion" and
#"missing" of the object are subsetted and updated appropriately, the offending genotypes
#are removed using the R/ASMap function subsetCross(). The linkage map can then be
#reconstructed.
pg <- profileGen(orderedmap, bychr = F, stat.type = c("xo", "dxo", "miss"),
                 id = "Line", xo.lambda = 14, layout = c(1, 3), lty = 2, cex = 0.7)
orderedmap.5 <- subsetCross(orderedmap, ind = !pg$xo.lambda)
orderedmap.6 <- mstmap(orderedmap.5, bychr = TRUE, dist.fun = "kosambi", trace = TRUE,
                       p.value = 1e-7, id = "Line")
chrlen(orderedmap.6)

profileMark(orderedmap.6, stat.type = c("seg.dist", "prop", "dxo", "recomb"),
            layout = c(1, 5), type = "l")

profileMark(orderedmap.6, stat.type = c("seg.dist", "dxo", "erf", "lod"), id = ("Line"),
            layout = c(1, 5), type = "l")

#pull the map and print it out
pull.map(orderedmap.6)
write.cross(cross=orderedmap.6, format = "csv", filestem = "andrea_RILs_map")
  ## note - exported the map to excel:
      #manually removed the linkage groups that were not chromosomes, renamed to andrea_RILs_map_justchromosomes.csv

##### QTL analysis for anthocyanin and flavonol concentration in RILs pop #####

setwd("~/RILs/")

library(qtl)

#LOAD THE DATA (Data available upon request, please contact cris.kuhlemeier@ips.unibe.ch)
DATA <- read.cross ("csv","~/RILs/", "axex_RILs_map_forpub.csv", crosstype = "riself")

DATA <- convert2riself(DATA)
# the default for a cross is backcross or intercross
# if it is a RIL set, we need to convert the DATA using convert2riself

#to visualize missing data and markers positions on a map
plotMissing(DATA)
plotMap(DATA)


jittermap(DATA) 
# some markers are often at the same position, use jittermap

DATA <- calc.genoprob (DATA) 
# calculate the probability that the observed genotype is the true genotype by taking into account a small genotyping error rate

geno.image (DATA)
# plot the colored genotype table

summaryMap(DATA)

#### quick plots ####

plot(DATA$pheno$anthocyanins.limb, DATA$pheno$flavonols.limb, xlab = c("Anthocyanin absorbance"), ylab = c("Flavonol absorbance"), ylim = c(0,100), xlim = c(0, 60))
hist(DATA$pheno$anthocyanins.limb)
hist(DATA$pheno$flavonols.limb)

#### PLOT QTL for anthos and flavonols, regular scanone ####
## anthos
#regular scanone
IMantho=scanone(DATA,pheno.col=3)
# to calculate the threshold, we perform 1000 permutation of the scanone (re-run with 10000 for better estimation, 1k is fine for quick check)
perm.antho=scanone(DATA,pheno.col=6,n.perm=1000)
summary(perm.antho)
#LOD thresholds (1000 permutations)
#    lod
#5%  2.60
#10% 2.27
# give the lod thresholds obtained at 5% and 10%
antho.lod=summary(perm.antho, alpha=0.05)
antho.lod=2.60
plot(IMantho) #red = black here, IM.fl = IM.fl.2part
abline(h = 2.60, lty=2)

#flavonols are bimodal, use 2part model
IM.fl.np <- scanone(DATA, pheno.col = 2, model = c("np")) #try nonparametric
IM.fl.2part <- scanone(DATA, pheno.col = 2, model = c("2part")) #try 2part model
IM.fl=scanone(DATA,pheno.col=2) #standard model
# to calculate the threshold, we perform 1000 permutation of the scanone (do 10k for good estimate)
perm.fl=scanone(DATA,pheno.col=2,n.perm=1000)
summary(perm.fl)
#LOD thresholds (1000 permutations)
#    lod
#5%  2.55
#10% 2.31
# give the lod thresholds obtained at 5% and 10%
flavonol.lod=summary(perm.fl, alpha=0.05)
flavonol.lod=2.55
plot(IM.fl, IM.fl.np, IM.fl.2part, col = c("red", "blue", "black")) #red = black here, IM.fl = IM.fl.2part
plot(IM.fl.np, IM.fl.2part, col = c("red", "blue", "black"))
#I think use 2part analysis since data is sort of bimodal

#plot together
plot(IMantho,IM.fl,col=c("red","black"), ylim=c(0,85))
legend("topright",c("anthocyanins","flavonols"),col=c("red","black"),lty=1,lwd=2)
#abline(h=0)
abline(h=2.6,lty=2)	

#plot together
plot(IMantho,IM.fl.2part,col=c("red","black"))
legend("topright",c("anthocyanins","flavonols"),col=c("red","black"),lty=1,lwd=2)
#abline(h=0)
abline(h=2.6,lty=2)	



#### PVE for antho and flavonol QTLs ####
#get LOD of each peak and sample size

#### PVE for anthos and flavonols ###
##anthos
chr2 <- IMantho[ which(IMantho$chr =='2'), ]
head(chr2[order(-chr2$lod),], n = 20)
max(chr2)
#25.53844
#LOD for chr2 qtl = 25.68480
#n=189 (from DATA total n, and running scanone dropping one)
#forumula 1 - 10^(-(2/n)*LOD) 
1-(10^(-(2/189)*25.68480))
#46% PVE, 0.4651855

#getting lod/PVE for smaller qtl, chr 1
chr1 <- IMantho[ which(IMantho$chr=='1'), ]
head(chr1[order(-chr1$lod),], n = 20)
max(chr1)
#LOD for chr1 qtl = 7.610192
#n=189 (from DATA total n, and running scanone)
#forumula 1 - 10^(-(2/n)*LOD) 
1-(10^(-(2/189)*7.610192))
#17% PVE, 0.1692528

#getting lod/PVE for smaller qtl, chr 7
chr7 <- IMantho[ which(IMantho$chr=='7'), ]
head(chr7[order(-chr7$lod),], n = 20)
max(chr7)
#LOD for chr7 qtl = 3.298620
#n=189 (from DATA total n, and running scanone)
#forumula 1 - 10^(-(2/n)*LOD) 
1-(10^(-(2/189)*3.298620))
#8% PVE, [1] 0.07722893

#getting lod/PVE for smaller qtl(s), chr 3
chr3 <- IMantho[ which(IMantho$chr=='3'), ]
head(chr3[order(-chr3$lod),], n = 20)
max(chr3)
#max LOD for chr3 qtl = 2.709738
#n=189 (from DATA total n, and running scanone)
#forumula 1 - 10^(-(2/n)*LOD) 
1-(10^(-(2/189)*2.709738))
#7% PVE, [1] 0.0685163

#PVE for flavonols for chr 2
#consider overall LOD score, which is LOD(π,μ) or lod.p.mu
chr2.fl <- IM.fl.2part[ which(IM.fl.2part$chr =='2'), ]
head(chr2.fl[order(-chr2.fl$lod.p.mu),], n = 20)
#LOD for chr2 qtl = 79.92402
#n=189 (from DATA total n, and running scanone)
#forumula 1 - 10^(-(2/n)*LOD) 
1-(10^(-(2/189)*79.92402))
#87% PVE, 0.8573594

#PVE for flavonols for chr 7
#consider overall LOD score, which is LOD(π,μ) or lod.p.mu
chr7.fl <- IM.fl.2part[ which(IM.fl.2part$chr =='7'), ]
head(chr7.fl[order(-chr7.fl$lod.p.mu),], n = 20)
#LOD for chr7 flav qtl = 4.564099
#n=189 (from DATA total n, and running scanone)
#forumula 1 - 10^(-(2/n)*LOD) 
1-(10^(-(2/189)*4.564099))
#10.5% PVE, 0.105248

#PVE for flavonols for chr 4
#consider overall LOD score, which is LOD(π,μ) or lod.p.mu
chr4.fl <- IM.fl.2part[ which(IM.fl.2part$chr =='4'), ]
head(chr4.fl[order(-chr4.fl$lod.p.mu),], n = 20)
#LOD for chr4 flav qtl = 3.028875
#n=189 (from DATA total n, and running scanone)
#forumula 1 - 10^(-(2/n)*LOD) 
1-(10^(-(2/189)*3.028875))
#7.11% PVE, 0.07114395


#######  PLOT QTL WITH ALLELIC DIRECTION
############# PLOT QTL WITH ALLELIC DIRECTION ####
####run effectscan to get the allelic direction (needs to be done on sim.geno(DATA) to ge all the genotype)
#for an intercross, it will show the addititve (a) effects of the QTLs
# and so a= -1, 0, 1, for the AA, AB and BB genotypes. 
#but we have RILs:
#In R/qtl, RIL genotypes are encoded like AA/BB, and the effects are like (ave_BB - ave_AA)/2.
#So positive values indicate B allele increases the phenotype and negative values indicate A allele increases the phenotype.
#Thus in my case, anything -1 is axillaris (AA) and anything +1 is exserta (BB)
#from https://www.rdocumentation.org/packages/qtl/versions/1.46-2/topics/effectscan
find.pheno(DATA,pheno = c("anthocyanins.limb")) #3
#is normal?
hist((DATA$pheno$anthocyanins.limb))
#no not really - div by 2?
#hist((DATA$pheno$anthocyanins.limb)/2)

#try /2
#DATA$pheno$anthocyanins.limb <- (DATA$pheno$anthocyanins.limb)/2

#effect scan of additive effects
effectscan(sim.geno(DATA),pheno.col = 3)
#this will plot the allelic direction but we just want to convert the scanone object with LODscore * (+/-)
#no need for a graph : draw = FALSE
allelic = effectscan(sim.geno(DATA),pheno.col = 3, draw = FALSE)
effectscan(sim.geno(DATA), pheno.col = 3, chr="2", mtick="triangle", get.se=TRUE) #showing confidence intervals on chr 2

#run scanone
IMallelic=scanone(DATA,pheno.col=3)
#dropping 1 indiv with missing phenotype, so n = 194

#let's now change the allelic direction of the LOD score
#by dividing the allelic effect by the absolute value of the allelic effect, we get 1 or -1
#then multiply the LOD score by that and you get a directed LOD score
IMallelic[,3] = IMallelic[,3] * allelic[,3]/abs(allelic[,3])


#calculate significance/LOD threshold

perm=scanone(DATA,pheno.col=3,n.perm=1000)
# to calculate the threshold, we perform 1000 permutation of the scanone
summary(perm)
#LOD thresholds (1000 permutations)
#    lod
#5%  2.60
#10% 2.35
# give the lod thresholds obtained at 5% and 10%
antho.lod=summary(perm, alpha=0.05)
antho.lod=2.60

plot(IMallelic,col=c("black","red","blue"), ylim=c(-20,25))
legend("topright",c("anthos"),col=c("black"),lty=1,lwd=2)
abline(h=0)
abline(h=antho.lod,lty=2)	
abline(h=-antho.lod,lty=2)

plot(IMallelic,chr=7,col=c("black","red","blue"), ylim=c(-20,25))
legend("topright",c("anthos"),col=c("black"),lty=1,lwd=2)
abline(h=0)
abline(h=antho.lod,lty=2)	
abline(h=-antho.lod,lty=2)


### plot allelic qtl with both anthos and flavonols####

find.pheno(DATA,pheno = c("flavonols.limb")) #2
#is normal?
hist((DATA$pheno$flavonols.limb))
#no, data very bimodal - use 2part algorithm since not normally distrib
hist(log(DATA$pheno$flavonols.limb))

#allelic flavonols, pheno col 2
IM.fl=scanone(DATA,pheno.col=2,model = c("2part")) #regular qtl scanone
IM.FL.allelic=scanone(DATA,pheno.col=2,model = c("2part"))
allelicFl = effectscan(sim.geno(DATA),pheno.col = 2, draw = FALSE)
IM.FL.allelic[,3] = IM.FL.allelic[,3] * allelicFl[,3]/abs(allelicFl[,3])

## plot together
#only possible to plot 3 scanone objects
plot(IMallelic,IM.FL.allelic,col=c("red","black"), ylim=c(-100,40))
legend("bottomright",c("anthocyanins","flavonols"),col=c("red","black"),lty=1,lwd=2)
abline(h=0)
abline(h=2.6,lty=2)	
abline(h=-2.6,lty=2)

plot(IMallelic,IM.FL.allelic,col=c("red","black"), ylim=c(-90,25))
legend("bottomright",c("anthocyanins","flavonols"),col=c("red","black"),lty=1,lwd=2)
abline(h=0)
abline(h=2.6,lty=2)	
abline(h=-2.6,lty=2)

plot(IMallelic,IM.FL.allelic,col=c("red","black"), ylim=c(-90,25))
#legend("bottomright",c("anthocyanins","flavonols"),col=c("red","black"),lty=1,lwd=2)
abline(h=0)
abline(h=2.6,lty=2)	
abline(h=-2.6,lty=2)


####everything above here is published in Berardi et al. ####