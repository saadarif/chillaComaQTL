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


#You could not continue with your QTL mapping approach as before but try with the cold6.clean cross
#Note I've put step to 0 so that it does not estimate any intermarker distances as this does not make much sense
#if you haven't estimated a genetic map. I would recommend that in additional to the analysis with no estimated genetic map
#you could do one with a genetic map estimated assuming F.gen=5, you can covert the current cross to F5 and simulatenously estimate the map as follows:
#cold6.clean.F5<- convert2bcsft(cold6.clean, F.gen = 5)
#Have a look at the genetic map before proceeding? Does it look reasonable?

#------------------------------------------------------------------------------------------------------------------------------------------
# Convert the data to F5 and estimate a genetic map

cold6.clean.F5<- convert2bcsft(cold6.clean, F.gen = 5)
plotMap(cold6.clean.F5, alternate.chrid = T)
summaryMap(cold6.clean.F5)
#         n.mar length ave.spacing max.spacing
# s12613     17    4.8         0.3         1.1
# s12902      2    0.5         0.5         0.5
# s12903     18   11.1         0.7         2.4
# s12911     43   64.8         1.5        11.7
# s12916    179   92.8         0.5        12.2
# s12929      9   26.0         3.2        16.2
# s12943      1    0.0          NA          NA
# s12947      3    0.9         0.4         0.9
# s13010      2    1.3         1.3         1.3
# s13047     32   17.7         0.6         9.0
# s13060      1    0.0          NA          NA
# s13077      2    3.5         3.5         3.5
# s13082     19    6.5         0.4         2.7
# s13088     11    5.2         0.5         1.6
# s13099      7    1.6         0.3         0.5
# s13117     63   95.2         1.5        40.8
# s13228      1    0.0          NA          NA
# s13230      1    0.0          NA          NA
# s13248     19    8.2         0.5         1.7
# s13250      4    2.5         0.8         1.7
# s13266    200  145.4         0.7        55.5
# s13334     39   19.9         0.5         1.8
# s13335      9   13.1         1.6         6.9
# s13337    288  256.7         0.9        19.8
# s13339      1    0.0          NA          NA
# s13340    419  179.6         0.4         6.0
# s13499      6    2.1         0.4         1.1
# s13767      3    2.3         1.2         1.3
# s3828       1    0.0          NA          NA
# overall  1400  962.0         0.7        55.5

est.rf(cold6c.F5)                                                       
lg<-formLinkageGroups(cold6c.F5, max.rf = 0.35, min.lod = 6) 
table(lg[,2])
#   1    2    3    4    5 
#1167  125   79   27    2 

checkAlleles(cold6.clean.F5)
# checkAlleles not available for cross type bcsft.

cold6 <- calc.genoprob(cold6, step = 0)   # calculate additional information: QTL genotype probabilities, step = density of the grin (in cM) 

######## scanone

result_em_cold6 <- scanone(cold6, pheno.col=1, method= "em")
result_hk_cold6 <- scanone(cold6, pheno.col=1, method= "hk")
result_ehk_cold6 <- scanone(cold6, pheno.col=1, method= "ehk")
### find significance threshold
perm_em_cold6<-scanone(cold6, method="em", n.perm = 1000)
perm_hk_cold6<-scanone(cold6, method="hk", n.perm = 1000)
perm_ehk_cold6<-scanone(cold6, method="ehk", n.perm = 1000)
summary(perm_em_cold6.clean.F5)
#lod
#5%  4.66
#10% 4.18
summary(perm_hk_cold6.clean.F5)
#lod
#5%  4.00
#10% 3.58
summary(perm_ehk_cold6.clean.F5)  
#lod
#5%  153
#10% 153
summary(result_em_cold6.clean.F5, threshold=3.5, df=TRUE)
#         chr     pos  lod
#22061 s13337  0.0839 6.81
#23306 s13340 80.2843 5.12

### Find QTL support intervals
lodint(result_em_cold6c.F5, chr="s13337")
#          chr      pos      lod
#22061  s13337 0.083871 6.811770
#22061  s13337 0.083871 6.811770
#20598  s13337 9.233870 5.119635
lodint(result_em_cold6c.F5, chr="s13340")
#         chr      pos      lod
#26993 s13340 27.80865 3.418549
#23306 s13340 80.28428 5.124599
#23639 s13340 96.77336 3.613436
# Which markers "flank" the intervals?
find.markerpos(cold6a, "22061")     
#         chr      pos
#22061 s13337 0.083871                
find.markerpos(cold6a, "20598")     
#         chr     pos
#20598 s13337 2.26785
find.markerpos(cold6a, "26993")       
#       chr      pos
#26993 s13340 5.707598
find.markerpos(cold6a, "23639")       
#         chr      pos
#23639 s13340 12.90732

######## scantwo

cold6.clean.F5 <- calc.genoprob(cold6.clean.F5, step = 0) 
result_scan2_hk_cold6.clean.F5 <- scantwo(cold6.clean.F5, pheno.col=1, method= "hk", clean.output=TRUE)       
perm_scan2_hk_cold6.clean.F5<- scantwopermhk(cold6.clean.F5, n.perm = 1000)  # takes a while
scan2_hk_cold6.clean.F5_summary <- summary(result_scan2_hk_cold6.clean.F5, perms = perm_scan2_hk_cold6.clean.F5, alpha = 0.2, pvalues=T)
scan2_hk_cold6.clean.F5_summary
#                pos1f pos2f lod.full pval lod.fv1 pval lod.int pval     pos1a pos2a lod.add pval lod.av1 pval
#cs13337:cs13340  6.97  29.4       13    0    6.15 0.42    1.61    1      6.97  30.1    11.4    0    4.54 0.02

write.table(scan2_hk_cold6.clean.F5_summary, "../../Analysis_in_R/version6_scan2.csv", row.names = F)


######## Stepwise MQM 

pen6c<-calc.penalties(perm_scan2_hk_cold6.clean.F5, alpha = 0.05)
pen6c
#    main    heavy    light 
#4.069643 6.644801 3.246213 
cold6.clean.F5<- sim.geno(cold6.clean.F5, step = 1, n.draws = 32, error.prob = 0.0001)
step_qtl_out_6c <- stepwiseqtl(cold6.clean.F5, verbose=T, max.qtl=5, penalties = pen6c, method = "imp")
summary(step_qtl_out_6c)
#    QTL object containing imputed genotypes, with 32 imputations. 
# 
# name    chr       pos n.gen
# Q1 s12916@16.7 s12916 16.747634     3
# Q2  s13337@0.1 s13337  0.083871     3
# Q3 s13340@30.1 s13340 30.141236     3
# 
# Formula: y ~ Q1 + Q2 + Q3 + Q1:Q3 
# 
# pLOD:  4.445 
out1_6c <- fitqtl(cold6.clean.F5, qtl = step_qtl_out_6c, method = "imp", formula = y ~ Q1 + Q2 + Q3 + Q1:Q3, get.ests = T )
summary_out1_6c<-summary(out1_6c)
summary_out1_6c

# fitqtl summary
# 
# Method: multiple imputation 
# Model:  normal phenotype
# Number of observations : 86 
# 
# Full model result
# ----------------------------------  
#   Model formula: y ~ Q1 + Q2 + Q3 + Q1:Q3 
# 
#       df        SS         MS      LOD     %var Pvalue(Chi2)    Pvalue(F)
# Model 10 12087.273 1208.72732 19.90036 65.54907 2.553513e-15 9.170442e-14
# Error 75  6352.764   84.70352                                            
# Total 85 18440.037                                                       
# 
# 
# Drop one QTL at a time ANOVA table: 
#   ----------------------------------  
#                           df Type III SS    LOD  %var F value Pvalue(Chi2) Pvalue(F)    
#   s12916@16.7              6        3783  8.724 20.51   7.443            0  2.85e-06 ***
#   s13337@0.1               2        4948 10.755 26.83  29.205            0  4.17e-10 ***
#   s13340@30.1              6        5876 12.229 31.86  11.561            0  4.01e-09 ***
#   s12916@16.7:s13340@30.1  4        2901  7.025 15.73   8.563            0  9.54e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# 
# Estimated effects:
#   -----------------
#                               est      SE      t
# Intercept                 48.0808  1.4441 33.294
# s12916@16.7a              -0.6019  1.2432 -0.484
# s12916@16.7d              -0.7837  2.7112 -0.289
# s13337@0.1a                9.1937  1.2640  7.274
# s13337@0.1d               -2.2274  2.3624 -0.943
# s13340@30.1a               1.7529  1.3843  1.266
# s13340@30.1d              -5.2418  2.7514 -1.905
# s12916@16.7a:s13340@30.1a  3.6557  1.4797  2.471
# s12916@16.7d:s13340@30.1a -6.2630  2.7269 -2.297
# s12916@16.7a:s13340@30.1d  2.5932  2.5565  1.014
# s12916@16.7d:s13340@30.1d 22.6301  5.4612  4.144


plot(step_qtl_out_6c, alternate.chrid = T)



