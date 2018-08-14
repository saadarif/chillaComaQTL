# QTL mapping of cold tolerance in D. ananassae from Bangkok
# skript modified from S. Arif ("CheckvariantcallDist.R")
# markers generated with Stacks v. 1.45 ("QTL_Mapping_pipeline4")

                                 #### Version 6 ####

                                 
                                 
#SNPs: second corrections applied; filtered for heterozygozity                                 
                                 
                                 
setwd("/media/koeniger/Elements1/RAD-Seq/QTL_Mapping_pipeline_4/stacks_out/cor_stacks_2")

#Step 1 Filtering on expected genotype and allele frequencies

#We will use markers that were called in 80% (75/94) of the progeny
snps6 <- read.table("batch_1.genotypes_75.onemap.txt", skip=3, na.strings=c("-","--"), stringsAsFactors = F)

#lowercase h is the original genotype call for heterozygote and captial H is the correction applied by the genotypes utility in stacks
#convert everything to uppercase
numsnps6= nrow(snps6)  #3092
numind6=ncol(snps6)-1  #94

library("stringr")
snps6[,2:ncol(snps6)] <- data.frame(lapply(snps6[,2:ncol(snps6)], function(v) {
  return(as.character(toupper(v)))
}))

head(snps6)

#make a dataframe to each marker name, heterozygity (total het call over all calls) and allele frequency drift (absoulte different in a-b allele fre)
marker_summ6 <- data.frame(markers=snps6$V1, het=numeric(length=numsnps6), drift=numeric(length=numsnps6), PropMiss=numeric(length=numsnps6),
                           numAA=numeric(length=numsnps6), numAB=numeric(length=numsnps6), numBB=numeric(length=numsnps6))


for (i in 1:numsnps6){
  #do the genotype counts
  x <- table(t(snps6[i,2:ncol(snps6)]), useNA = "always")
  #calculate heterozygote frequencies
  marker_summ6$het[i] <- x[3]/sum(x[1:3])
  #calculate the allele frequency drift
  #1. calculate allele frequency of one allele
  a <- ((2*x[1])+x[3])/(2*sum(x[1:3]))
  #2. then the absolute differenece of it's frequency from 0.5
  marker_summ6$drift[i] <- abs(0.5-a)
  #finally addinformation on proportion of missing indivduals for the tag
  marker_summ6$PropMiss[i] <- x[4]/numind6
  marker_summ6$numAA[i] <- x[1]
  marker_summ6$numAB[i] <- x[3]
  marker_summ6$numBB[i] <- x[2]
}

head(marker_summ6)

# noow use the criter for het >= 0.15 and <= 0.35 and drift <=0.1 to eliminate crazy marker
# experiment with marker het between 0.1-0.4 and drift upto 0.2
marker_clean6 <- marker_summ6[which((marker_summ6$het >= 0.15 & marker_summ6$het <= 0.35) & marker_summ6$drift <= 0.1),]
nrow(marker_clean6) # 1400
head(marker_clean6)
#Alternatively using segregation distortion tests
#The following can be use Instead to check against a segregation distortion of 0.375:0.25:0.375
#Calculate all the p-values first
#marker_seg6 <- data.frame(markers=snps6$V1, pvals=numeric(length=numsnps6), adjPvals=numeric(length=numsnps6))
#exp_p <- c(0.375,0.25,0.375)
#for (i in 1:numsnps6){
  #avoid missing genotypes for markers
#  if(any(is.na(c(marker_summ6$numAA[i],marker_summ6$numAB[i],marker_summ6$numBB[i])))) {
#    marker_seg6$pvals[i] = NA
#  }
#  else{
    #print(c(marker_summ$numAA[i],marker_summ$numAB[i],marker_summ$numBB[i]))
#    test_results6 <-chisq.test(c(marker_summ6$numAA[i],marker_summ6$numAB[i],marker_summ6$numBB[i]), p=exp_p  )
#    marker_seg3$pvals[i] <- test_results3$p.value
#  }
#}

#Adjust the p-values usign Bejnaminin-Hochberg, could use bonferroni to be more conservative
#marker_seg3$adjPvals <- p.adjust(marker_seg3$pvals, method="BH")

#Filter markers based on NAs and padj <0.05  
#marker_clean_seg3 <- marker_summ6[which(marker_seg3$adjPvals > 0.05),]
#nrow(marker_clean_seg3)


# ----------------------------------------------------------------------------------------------------------

#Step 2: Prepare an input file for R/QTL from the input file

#we'll need the phenotypes
pheno <- read.csv("../../Analysis_in_R/pheno_file.csv")

#construct the first column, the pheno colun, for the read.cross file
pheno_cross <- c("chill_coma", " ", " ", pheno$pheno)

#now we need scaffold and positions information, this is stored "batch_1.catalog.tags.tsv.gz"
# output file from the stacks genotypes run
#system("gunzip batch_1.catalog.tags.tsv.gz")

mark_info6 <- read.table("batch_1.catalog.tags.tsv", header=F)
head(mark_info6)
str(mark_info6)
#catalog id from the one map file is in the 3rd column (it is the row number), we want to recover the 4th and the 5th col (scaffold position in kb)
#if filtered above
cat_ids6 <- as.numeric(str_sub(marker_clean6$markers, 2, -1))
#if NOTfiltered
#cat_ids6 <- as.numeric(str_sub(snps6$V1, 2, -1))

marker_header6 <- mark_info6[cat_ids6, 3:5]
str(marker_header6)
#convert all genotypes calls to same case - already done on lines 24-26
#recoded <- apply(snps[2:ncol(snps)],2,tolower)

#replace the marker columns in the snps file with the marker name, scaffold and position info
#if filtered
new_snps6 <- cbind(marker_header6, snps6[which(snps6$V1 %in% marker_clean6$markers), 2:ncol(snps6)])
#if NOT filtered
#new_snps6 <- cbind(marker_header6, snps6[, 2:ncol(snps6)])
head(new_snps6)
#reorder the columns by scaffold than positions
new_snps6 <- new_snps6[with(new_snps6, order(V4, V5)), ]
#transpose the new_snps and combine with the phen0o_cross
#save the file as .csv
crossfile6 <- cbind(pheno_cross, t(new_snps6))
str(crossfile6)
#getrid of the rownames and colnames
rownames(crossfile6)=NULL
colnames(crossfile6)=NULL
#save the file
write.table(crossfile6, "../../Analysis_in_R/Cold_r75_crossfile6.csv", quote=F, sep=",",col.names=F, row.names =F)
#Replace the long "scaffolds_" with "s" in the file 
system("sed -i 's/scaffold_/s/g' ../../Analysis_in_R/Cold_r75_crossfile6.csv")

#-------------------------------------------------------------------------------------------------------
#  Step 3: Playing around in R/Qtl

#Note this assumes the cross file is in the same directory and the working directory is the the one where the script is stored
#change accordingly

cold6 <- read.cross(format="csv", ".", "Cold_r75_crossfile6.csv", genotypes=c("A",  "H", "B"), estimate.map = F)
#I have not specified F.Gen as i will not estimating a genetic map
#you may want to repeat the analysis with cleaned version below and genetic distances estimated assuming F.gen=5, for comparison
#you can convert your cross object to F.gen =t with estimated map as follows

#cold6<- convert2bcsft(cold6, F.gen = 5)

summary(cold6)
#rescale map from kb to mb
cold6 <- rescalemap(cold6, 1e-6) 

#----------------------------------------------------------------------------------------------------------------------------
#Check genotyping quality after filtering
geno.image(cold6)
#in an ideal world this should only be contiguous blocks of red,green (homozygous) and some blue (heterozygous)

#We still need to clean this up a bit
#R/qtl does deal with genotyping error but as should be obvious from the plot above
#there is still quite a bit of error here.

#We could do two things:
#(i)get rid of problematic indivdiuals
#(ii) get get rid od problematic markers in singl individuals.

#In this cleaning approach both filtering methods rely on detecting unexpected double crossover events

#The following is adapted from pages 33-39 from this tutorial
#http://www.rqtl.org/tutorials/geneticmaps.pdf

#Lets remove problematic individuals first as they will make detecting genotyping errors more tenable later

#Plot the number of double crossovers per individual
plot(countXO(cold6), ylab="Number of crossovers")
#most lines have double crossovers betweem 30-150 but ehre are clear outliers. Lets remove these for now
cold6 <- subset(cold6, ind=(countXO(cold6) < 200))
#We lose 8 lines
nind(cold6)
#Did this qualitatively improve our genotyping?
geno.image(cold6)
#I'm not sure if retaining the indvidual with missing data is helpful or not, it is not helpful to have the individual for the nexy
#step (identifying potential genotyping errors ) so I will remove that individual
cold6 <- subset(cold6, ind=(ntyped(cold6)>700))

#Now can We identify potentially erroneous double crossovers assumign a genotyping error rate of 1%
cold6 <- calc.errorlod(cold6, error.prob=0.01)

#We will remove any marker with an abritrary large LOD score, for e.g. greater than 5 (cutoff)
nrow(toperr <- top.errorlod(cold6, cutoff=5))

#remove markers with LOD scores greater than five and store it as a new cross object
cold6.clean <-cold6
for(i in 1:nrow(toperr)) {
  chr <- toperr$chr[i]
  id <- toperr$id[i]
  mar <- toperr$marker[i]
  cold6.clean$geno[[chr]]$data[id, mar] <- NA
}

#NOTE: I would recommend looking at the help page for calc.errorlod to find out what is happening and also looking at the function
#clean.geno() as an alternative/complimentary

summary(cold6)
summary(cold6.clean)
#We haven't lost a lot of markers but we have lost 9 indiviuals

#Has the genotyping qualtiatively improved?
geno.image(cold6.clean)



#------------------------------------------------------------------------------------------------------------------------------------------
#YOU COULD IGNORE THE FOLLOWING IF YOU LIKE, it is not directly related to above it is just checking the genetic map
#estimating the map using the est.map function
newmap <- est.map(cold6, n.cluster=6)
plotMap(cold6, newmap, alternate.chrid=T)
#There seem to be some potentially problematic markers on s13337  (and maybe the terminal one on s13117)
#Lets calculate pairwise recombination frequency
#Can we remove the problematic (unlinked markers) on s13337
#cold6 <- est.rf(cold6)
#and plot it 
#plotRF(cold6, alternate.chrid=T) 
#checkAlleles(cold6)
#cold6<-drop.markers(cold6, c("19153", "19632"))
#cold6 <- est.rf(cold6)
#plotRF(cold6, alternate.chrid=T) 
#PLEASE STOP IGNORING NOW
#----------------------------------------------------------------------------------------------------------------------------------------------

#You could not continue with your QTL mapping approach as before but try with the cold6.clean cross
#Note I've put step to 0 so that it does not estimate any intermarker distances as this does not make much sense
#if you haven't estimated a genetic map. I would recommend that in additional to the analysis with no estimated genetic map
#you could do one with a genetic map estimated assuming F.gen=5, you can covert the current cross to F5 and simulatenously estimate the map as follows:
#cold6.clean.F5<- convert2bcsft(cold6.clean, F.gen = 5)
#Have a look at the genetic map before proceeding? Does it look reasonable?


cold6 <- calc.genoprob(cold6, step = 0)   # calculate additional information: QTL genotype probabilities, step = density of the grin (in cM) 
result_em_cold6 <- scanone(cold6, pheno.col=1, method= "em")
result_hk_cold6 <- scanone(cold6, pheno.col=1, method= "hk")
result_ehk_cold6 <- scanone(cold6, pheno.col=1, method= "ehk")
# find significance threshold
perm_em_cold6<-scanone(cold6, method="em", n.perm = 1000)
perm_hk_cold6<-scanone(cold6, method="hk", n.perm = 1000)
perm_ehk_cold6<-scanone(cold6, method="ehk", n.perm = 1000)
summary(perm_em_cold6)
#     lod
#5%  3.90
#10% 3.52
summary(perm_hk_cold6)
#     lod
#5%  3.94
#10% 3.53
summary(perm_ehk_cold6)  
# 167
summary(result_em_cold6, threshold=3.5, df=TRUE)
#                chr    pos  lod
#377          s12903  0.295 3.54   # not sig at 5%
#cs13337.loc2 s13337  2.084 7.39
#23306        s13340 11.450 5.78
summary(result_hk_cold6, threshold=3.5, df=TRUE)
#                chr    pos  lod
#  377        s12903  0.295 3.61
#cs13337.loc2 s13337  2.084 7.17
#23306        s13340 11.450 5.86

lodint(result_em_cold6, chr="s13337") 
#                chr      pos      lod
#22061        s13337 0.083871 6.896617
#cs13337.loc2 s13337 2.083871 7.394478
#21689        s13337 6.757392 5.248228
lodint(result_em_cold6, chr="s13340") 
#                 chr      pos      lod
#cs13340.loc10 s13340 10.05311 4.162521
#23306         s13340 11.45033 5.777952
#23858         s13340 13.69258 3.543581

######## scantwo

cold6 <- calc.genoprob(cold6, step = 0) 
result_scan2_hk_cold6 <- scantwo(cold6, pheno.col=1, method= "hk", clean.output=TRUE)       
summary(result_scan2_hk_cold6)
perm_scan2_hk_cold6<- scantwopermhk(cold6, n.perm = 1000)  # takes a while
scan2_hk_cold6_summary <- summary(result_scan2_hk_cold6, perms = perm_scan2_hk_cold6, alpha = 0.2, pvalues=T)
#                pos1f pos2f lod.full  pval lod.fv1  pval lod.int pval     pos1a pos2a lod.add pval lod.av1 pval
#cs13337:cs13340  2.08  13.1     12.2 0.001    5.07 0.951    1.03    1      2.08  10.1    11.2    0    4.04 0.03

write.table(scan2_hk_cold6_summary, "../../Analysis_in_R/version6_scan2.csv", row.names = F)


######## MQM by imputation

#Multiple QTL model:

qtl_model6 <- sim.geno(cold6, step = 1, n.draws = 32, error.prob = 0.0001)
qtl6 <- makeqtl(qtl_model6,chr=c("s13337", "s13340"), pos=c(2.084, 11.450)) # candidate regions
qtl6

# QTL object containing imputed genotypes, with 32 imputations. 
# 
#         name    chr     pos n.gen
#Q1  s13337@2.1 s13337  2.0839     3
#Q2 s13340@11.5 s13340 11.4503     3

plot(qtl6, alternate.chrid=T, main= "Genetic map version 6") # Map of the chromosomes with candidate loci

# possible model without interaction:

out.fqa6 <- fitqtl(qtl_model6, qtl=qtl6, get.ests=T, pheno.col=1,formula=y~Q1+Q2)
mqm_imp_cold6_summary <- summary(out.fqa6)

# fitqtl summary
# 
# Method: multiple imputation 
# Model:  normal phenotype
# Number of observations : 94 
# 
# Full model result
# ----------------------------------  
#   Model formula: y ~ Q1 + Q2 
# 
#       df        SS        MS      LOD     %var Pvalue(Chi2)    Pvalue(F)
# Model  4  8122.773 2030.6932 11.19068 42.20365 1.725543e-10 5.021195e-10
# Error 89 11123.838  124.9869                                            
# Total 93 19246.611                                                      
# 
# 
# Drop one QTL at a time ANOVA table: 
#   ----------------------------------  
#             df Type III SS   LOD  %var F value Pvalue(Chi2) Pvalue(F)    
# s13337@2.1   2        3450 5.515 17.93   13.80            0  6.01e-06 ***
# s13340@11.5  2        2558 4.225 13.29   10.23            0     1e-04 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# 
# Estimated effects:
#   -----------------
#                  est      SE      t
# Intercept    48.5213  1.3745 35.302
# s13337@2.1a   7.5859  1.4518  5.225
# s13337@2.1d  -0.6892  2.6374 -0.261
# s13340@11.5a  5.9225  1.4428  4.105
# s13340@11.5d -4.7136  2.6821 -1.757

mqm_imp_cold6_tab <- mqm_imp_cold6_summary$result.drop 
write.table(mqm_imp_cold6_tab, file="mqm_imp_cold6.csv", sep="\t", quote=FALSE, row.names=TRUE) 

# alternatively: model with interaction:

out.fqa6a <- fitqtl(qtl_model6, qtl=qtl6, get.ests=T, pheno.col=1,formula=y~Q1*Q2)
mqm_imp_cold6a_summary <- summary(out.fqa6a)


# fitqtl summary
# 
# Method: multiple imputation 
# Model:  normal phenotype
# Number of observations : 94 
# 
# Full model result
# ----------------------------------  
#   Model formula: y ~ Q1 + Q2 + Q1:Q2 
# 
# df        SS        MS      LOD     %var Pvalue(Chi2)    Pvalue(F)
# Model  8  8509.512 1063.6890 11.91296 44.21304 4.697757e-09 2.340231e-08
# Error 85 10737.099  126.3188                                            
# Total 93 19246.611                                                      
# 
# 
# Drop one QTL at a time ANOVA table: 
#   ----------------------------------  
#   df Type III SS    LOD   %var F value Pvalue(Chi2) Pvalue(F)    
# s13337@2.1              6      3837.1 6.2368 19.937  5.0628        0.000  0.000175 ***
#   s13340@11.5             6      2944.6 4.9468 15.299  3.8851        0.001  0.001781 ** 
#   s13337@2.1:s13340@11.5  4       386.7 0.7223  2.009  0.7654        0.505  0.550683    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# 
# Estimated effects:
#   -----------------
#   est      SE      t
# Intercept                48.2083  1.4595 33.030
# s13337@2.1a               6.9167  1.7484  3.956
# s13337@2.1d               0.5324  2.9411  0.181
# s13340@11.5a              5.2267  1.6948  3.084
# s13340@11.5d             -4.9973  2.9084 -1.718
# s13337@2.1a:s13340@11.5a  2.2234  1.7519  1.269
# s13337@2.1d:s13340@11.5a -1.5691  3.5239 -0.445
# s13337@2.1a:s13340@11.5d -2.8247  3.5809 -0.789
# s13337@2.1d:s13340@11.5d  0.9469  6.0358  0.157

