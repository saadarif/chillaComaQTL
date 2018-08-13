#Looking at heterozygous calls in Annabella's data
#SA Feb 2018
#The data was out as onemap format from the genotypes utility in stacks
#This script is meant to flag markers that look suspicisous
#Markers/tags are eliminated based on the following criteria:
# 1. expected heterozgosity after 5 bouts of full sib mating is ~ 25% (Rumball elt a. 1994; Pollard 2012 (in Methods in Molecular biology))
# therefore we remove any markers with < 15% or > 35 % heterozygous calls
#2. Expected allele frequency drift (expected frequency of a and b alleles should be roughly 0.5) is ~ 0.03 for simulations closest to Annabella's line 

setwd("/media/data_disk/PROJECTS/Saad/Annabella/stacks.ref.annabella/cor_stacks_2/")

#-------------------------------------------------------------------------------------------------------------------------
#Step 1 Filtering on expected genotype and allele frequencies

#We will use markers that were called in 80% (75/94) of the progeny
snps <- read.table("batch_1.genotypes_75.onemap.txt", skip=3, na.strings=c("-","--"), stringsAsFactors = F)

#lowercase h is the original genotype call for heterozygote and captial H is the correction applied by the genotypes utility in stacks
#convert everything to uppercase
numsnps= nrow(snps)
numind=ncol(snps)-1

library("stringr")
snps[,2:ncol(snps)] <- data.frame(lapply(snps[,2:ncol(snps)], function(v) {
   return(as.character(toupper(v)))
}))

#make a dataframe to each marker name, heterozygity (total het call over all calls) and allele frequency drift (absoulte different in a-b allele fre)
marker_summ <- data.frame(markers=snps$V1, het=numeric(length=numsnps), drift=numeric(length=numsnps), PropMiss=numeric(length=numsnps),
                          numAA=numeric(length=numsnps), numAB=numeric(length=numsnps), numBB=numeric(length=numsnps))

for (i in 1:numsnps){
  #do the genotype counts
  x <- table(t(snps[i,2:ncol(snps)]), useNA = "always")
  #calculate heterozygote frequencies
  marker_summ$het[i] <- x[3]/sum(x[1:3])
  #calculate the allele frequency drift
  #1. calculate allele frequency of one allele
  a <- ((2*x[1])+x[3])/(2*sum(x[1:3]))
  #2. then the absolute differenece of it's frequency from 0.5
  marker_summ$drift[i] <- abs(0.5-a)
  #finally addinformation on proportion of missing indivduals for the tag
  marker_summ$PropMiss[i] <- x[4]/numind
  marker_summ$numAA[i] <- x[1]
  marker_summ$numAB[i] <- x[3]
  marker_summ$numBB[i] <- x[2]
}

# noow use the criter for het >= 0.15 and <= 0.35 and drift <=0.1 to eliminate crazy marker
#experiment with marker het between 0.1-0.4 and drift upto 0.2
marker_clean <- marker_summ[which((marker_summ$het >= 0.15 & marker_summ$het <= 0.35) & marker_summ$drift <= 0.1),]

#Alternatively using segregation distortion tests
#The following can be use Instead to check against a segregation distortion of 0.375:0.25:0.375
#Calculate all the p-values first
marker_seg <- data.frame(markers=snps$V1, pvals=numeric(length=numsnps), adjPvals=numeric(length=numsnps))
exp_p <- c(0.375,0.25,0.375)
for (i in 1:numsnps){
  #avoid missing genotypes for markers
  if(any(is.na(c(marker_summ$numAA[i],marker_summ$numAB[i],marker_summ$numBB[i])))) {
    marker_seg$pvals[i] = NA
  }
  else{
  #print(c(marker_summ$numAA[i],marker_summ$numAB[i],marker_summ$numBB[i]))
  test_results <-chisq.test(c(marker_summ$numAA[i],marker_summ$numAB[i],marker_summ$numBB[i]), p=exp_p  )
  marker_seg$pvals[i] <- test_results$p.value
  }
}

#Adjust the p-values usign Bejnaminin-Hochberg, could use bonferroni to be more conservative
marker_seg$adjPvals <- p.adjust(marker_seg$pvals, method="BH")

#Filter markers based on NAs and padj <0.05  
marker_clean_seg <- marker_summ[which(marker_seg$adjPvals > 0.05),]
#------------------------------------------------------------------------------------------------------------------------
#Step 2: Prepare an input file for R/QTL from the input file

#we'll need the phenotypes
pheno <- read.csv("/media/data_disk/PROJECTS/Saad/Annabella/Stacks_out_Rqtl_in/pheno_file.csv")

#construct the first column, the pheno colun, for the read.cross file
pheno_cross <- c("chill_coma", " ", " ", pheno$pheno)

#now we need scaffold and positions information, this is stored "batch_1.catalog.tags.tsv.gz"
system("gunzip batch_1.catalog.tags.tsv.gz")

#read in the catalog
mark_info <- read.table("batch_1.catalog.tags.tsv", header=F)
head(mark_info)
#catalog id from the one map file is in the 3rd column (it is the row number), we want to revover the 4th and the 5th col (scaffold position in kb)
#if filtered above
cat_ids <- as.numeric(str_sub(marker_clean_seg$markers, 2, -1))
#if NOTfiltered
#cat_ids <- as.numeric(str_sub(snps$V1, 2, -1))

marker_header <- mark_info[cat_ids, 3:5]

#convert all genotypes calls to same case - already done on lines 24-26
#recoded <- apply(snps[2:ncol(snps)],2,tolower)

#replace the marker columns in the snps file with the marker name, scaffold and position info
#if filtered
new_snps <- cbind(marker_header, snps[which(snps$V1 %in% marker_clean_seg$markers), 2:ncol(snps)])
#if NOT filtered
#new_snps <- cbind(marker_header, snps[, 2:ncol(snps)])

#reorder the columns by scaffold than positions
new_snps <- new_snps[with(new_snps, order(V4, V5)), ]
#transpose the new_snps and combine with the phen0o_cross
#save the file as .csv
crossfile <- cbind(pheno_cross, t(new_snps))

#getrid of the rownames and colnames
rownames(crossfile)=NULL
colnames(crossfile)=NULL
#save the file
write.table(crossfile, "Cold_r75_filtered2.csv", quote=F, sep=",",col.names=F, row.names =F)
#Replace the long "scaffolds_" with "s" in the file 
system("sed -i 's/scaffold_/s/g' Cold_r75_filtered2.csv")
#-------------------------------------------------------------------------------------------------------
#Step 3: Playing around in R/Qtl
library(qtl)
test <- read.cross(format="csv", ".", "Cold_r75_filtered2.csv", genotypes=c("A",  "H", "B"), estimate.map = F, F.gen = 5)

#rescale map from kb to mb
test <- rescalemap(test, 1e-6) 

#estimating the map using the est.map function
newmap <- est.map(test, error.prob = 0.01, n.cluster=4)
plotMap(test, newmap)

#plot of missing data
geno.image(test)

#one individual seesm to have lots of missign data, remove this individuaks
plot(test)
#number of typed markers per individual
ntyped(test)

#get rid of individuals with lots of missing markers
test <- subset(test, ind=(ntyped(test)>1100))

#duplicate individuals
cg <- comparegeno(test)
hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes")
rug(cg[lower.tri(cg)])
#look grrreat

#The chill comma phenotype doesn't look normal, consider transformation?


#Look at the reombination frequencies
test <- est.rf(test)
plot.rf(test, alternate.chrid = T)
checkAlleles(test, threshold=5)
#Note the error terms, some markers may need switching for be dropped

#an Example of multiple qtl mapping with augmenting missing data
#Use augmentation to fill the missing genotypes
test_aug <- mqmaugment(test, minprob = 0.1, verbose=T)
#alternatively impyting genos
test_fill<- fill.geno(test)


geno.image(test_aug)
geno.image(test_fill)
nind(test_aug)   #The augmented set contains more individuals then
nind(test)  #The unaugmented dataset
#So as a rule of thumb we can place (# individuals - 20) cofactors, which will be used in backward elimination.
cofactors <- mqmautocofactors(test_aug, 73, plot=T)
qtlresult <- mqmscan(test_aug, cofactors=cofactors, verbose=T, plot=T)
#Weird results

#try scanone
test <- calc.genoprob(test)
result_em <- scanone(test, pheno.col=1, method= "em")
#find highest peak
max(result_em)
#and the closest marker
find.marker(test, "scaffold_13337", 1.34)
#pull genotypes for this marker and use as covairate
g <- pull.geno(test)[,"10552"]
#impute the missing data and run as covariate



result_em_2 <- scantwo(test, pheno.col=1, method= "em")


#convert the data to a ril
test_ril <- convert2risib(test)
summary(test_ril)

test_ril<- calc.genoprob(test_ril)
result_em_ril <- scanone(test_ril, pheno.col=1, method= "em")

#stepwise qtl
#additive model only
step_qtl2 <- stepwiseqtl(test_fill, verbose=T, max.qtl=6, additive.only = T)
#with interaction
step_qtl_i2 <- stepwiseqtl(test, verbose=T, max.qtl=6)
#--------------------------------------------------------------------------------------------------

