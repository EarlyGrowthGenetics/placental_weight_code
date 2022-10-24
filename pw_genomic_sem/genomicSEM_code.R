###################
### GENOMIC SEM ###
###################


## Installing and loading packages: 
install.packages("devtools")
install_github("MichelNivard/GenomicSEM")
library(devtools)
library(GenomicSEM)


#STEP 1: MUNGE THE FILES 
files<-c("Clean_pw_fetal_sex_gest",
         "Clean_pw_maternal_sex_gest",
         "Clean_pw_paternal_sex_gest")

hm3<-"w_hm3.noMHC.snplist"  	## File provided by GenomicSEM developers
trait.names<-c("FetalPlacenta_sex_gest","MaternalPlacenta_sex_gest","PaternalPlacenta_sex_gest") 	## Name of traits in the model
N<-c(65405,61228,52392) 	## Sample size of traits in model

munge(files=files,hm3=hm3, trait.names=trait.names, N=N)

#STEP 2: RUN LD-SCORE REGRESSION 
traits<- c("FetalPlacenta_sex_gest.sumstats.gz", "MaternalPlacenta_sex_gest.sumstats.gz","PaternalPlacenta_sex_gest.sumstats.gz")
sample.prev <- c(NA,NA,NA)  	## Set to NA when continuous trait 
population.prev <- c(NA,NA,NA) 	## Set to NA when continuous trait 
 
ld <- "eur_w_ld_chr/" 		## File provided by GenomicSEM developers
wld <- "eur_w_ld_chr/" 		## File provided by GenomicSEM developers
trait.names<-c("FetalPlacenta_sex_gest","MaternalPlacenta_sex_gest","PaternalPlacenta_sex_gest")

PW_TRIOS<-ldsc(traits=traits, sample.prev=sample.prev, population.prev=population.prev, ld=ld, wld=wld,trait.names=trait.names)

save(PW_TRIOS, file="PW_TRIOS_sex_gest.RData")

## STEP 3: SPECIFY MODEL

covstruc<-PW_TRIOS

#Model:
Placenta <- '
Fetal =~ .5*MaternalPlacenta_sex_gest 	## make a latent fetal variable that loads on Maternal (with .5)
Fetal =~ .5*PaternalPlacenta_sex_gest 	## make a latent fetal variable that loads on Paternal (with .5)
Fetal =~ 1*FetalPlacenta_sex_gest 	## make a latent fetal variable that loads on Fetal (with 1)

Maternal =~ 1*MaternalPlacenta_sex_gest ## make a latent maternal variable that loads on Maternal with 1
Maternal =~ .5*FetalPlacenta_sex_gest	## make a latent maternal variable that loads on Fetal with.5

Paternal =~ 1*PaternalPlacenta_sex_gest	## make a latent paternal variable that loads on Paternal with 1
Paternal =~ .5*FetalPlacenta_sex_gest	## make a latent paternal variable that loads on Fetal with.5

MaternalPlacenta_sex_gest ~~ 0*MaternalPlacenta_sex_gest + 0*FetalPlacenta_sex_gest + 0*PaternalPlacenta_sex_gest ##no residual (co)variance between summary statistics
PaternalPlacenta_sex_gest ~~ 0*FetalPlacenta_sex_gest + 0*PaternalPlacenta_sex_gest	##no residual (co)variance between summary statistics
FetalPlacenta_sex_gest ~~ 0*FetalPlacenta_sex_gest	##no residual (co)variance between summary statistics

Fetal ~~ cov12*Maternal 	## free variances between latent variables
Fetal ~~ cov13*Paternal		## free variances between latent variables
Maternal ~~ cov23*Paternal	## free variances between latent variables

Fetal ~~var1*Fetal		## free variances of the latent variables
Maternal ~~ var2*Maternal	## free variances of the latent variables
Paternal ~~ var3*Paternal	## free variances of the latent variables

corFetal_Maternal := cov12/(sqrt(var1)*sqrt(var2)) 	## Calulated correlation between Fetal and Maternal latent variable
corFetal_Paternal := cov13/(sqrt(var1)*sqrt(var3))	## Calulated correlation between Fetal and Paternal latent variable
corMaternal_Paternal := cov23/(sqrt(var2)*sqrt(var3))	## Calulated correlation between Paternal and Maternal latent variable
'

estimation<-"DWLS"
std.lv=FALSE

results <- usermodel(covstruc=covstruc, model=Placenta,estimation=estimation,std.lv=std.lv)  	## Running genomicSEM model
results	## Printing results 