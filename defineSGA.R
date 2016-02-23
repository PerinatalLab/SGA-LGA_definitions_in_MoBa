####################################################################################
##### CALCULATE BIRTH-WEIGHT PERCENTILES (SKJAERVEN, HADLOCK, GARDOSI, MARSAL) #####

## This script uses four different methods to evaluate birth-weight.
## The output contains per-child info (gestational age, birthweight, sex), pregnancy and child IDs,
## and 8 new columns: percentiles and SGA indications (1/0) for each model.
## References for the model equations are found next to the formulas.


##### INPUT REQUIREMENTS #####

## This script requires three input files:
##
## 1. imputed maternal data, in the following format:
##      (INCLUDING HEADER LINE)
##      PREG_ID     AA85    AA86    AA87    flAA85  flAA86  flAA87
##      (this format is produced by script MoBa_imputing_MaternalHeightWeight.R)
##      AA85-87 can be taken directly from MoBa Q1
##      flAA85-87 are imputation flags. these flags are not required -
##      if you did not use the imputation script, you can turn it off in settings below
##
## 2. additional information directly from MFR, in the following format:
##      (NO HEADER LINE)
##      PREG_ID     CHILD_NUMBER    MOTHERS_AGE     GA(days)      SEX(1/2)     BIRTHWEIGHT(g)     PARITY
##      these columns correspond to the following variables from MFR:
##      PREG_ID_540     BARN_NR    MORS_ALDER     SVLEN_DG      KJONN     VEKT     PARITET_5
##
## 3. population means and standard deviations for Skjaerven's method, in the following format:
##      (INCLUDING HEADER LINE)
##      GA(completed weeks)       SEX(1/2)       MEAN      SD
##      this file can be retrieved from PerinatalLab GitHub:
##      https://raw.githubusercontent.com/PerinatalLab/SGA-LGA_definitions_in_MoBa/master/NSGA_BirthWeight_by_Skjaerven_2000_MEANSD.txt


##### SETTINGS #####

## PLEASE CHECK AND UPDATE THIS SECTION BEFORE PROCEEDING!

## do you have GitHub installed?
haveGit = TRUE

## if yes, enter its folder; if not, leave any value
gitDir = "~/Documents/gitrep/SGA-LGA_definitions_in_MoBa/"

## enter the folder which contains input data files (output will be written there as well)
file_dir = "~/Desktop/MoBa_v6/"

## the output file name will start with this:
outFileStem = "MOBA_PDB540_Linda_fetalWEIGHTandSGA"

### input file names:
## 1. this file contains the cleaned and imputed mother height/weight info
infile = "MOBA_PDB540_Linda_IMPUTED_maternalHgh1Wgh1Wgh2_20160223_32391fa.txt"
## does this file contain imputation flag columns?
imputeFlags = TRUE

## 2. this file contains the MFR data
mfrfile = "output_mfr_basicfetalinfo.csv"

## 3. this file contains the population data for Skjaerven's method (enter the complete path!)
skjaervenfile = "NSGA_BirthWeight_by_Skjaerven_2000_MEANSD.txt"


####################################################################################
####################################################################################

##### CALCULATIONS #####
## no changes should be needed here - just run everything below.

## install the SQLDF package if you do not have it yet
if(!require(dplyr)) install.packages('dplyr')
library(dplyr)
if(!require(gplots)) install.packages('gplots')
library(gplots)

## retrieve git hash if possible
if(haveGit){
        setwd(gitDir)
        hash = system(paste("git log --pretty=format:'%h' -n 1"),intern=TRUE)
} else { hash = "0000000" }

## get the date stamp
date_stamp = paste(unlist(strsplit(substr(Sys.time(),1,10),"-")), collapse="")
file_out = paste(file_dir, outFileStem, "_", date_stamp, "_", hash, ".txt",sep="")

## read input files
q1 = read.csv(paste(file_dir, infile, sep=""), sep="\t", header=T)
if(imputeFlags){
        names(q1) = c("PREG_ID","AA85","AA86","AA87","flAA85","flAA86","flAA87")
} else {
        names(q1) = c("PREG_ID","AA85","AA86","AA87")
}
head(q1); dim(q1)

mfr = read.csv(paste(file_dir, mfrfile, sep=""),sep=",", header=F)
names(mfr) = c("PREG_ID","CHILDNUM","AGE","GA","SEX","BIRTHWEIGHT","PARITY")
head(mfr); dim(mfr)

## merge them
### note that a unique fetal ID is also created - it can be useful to keep track of the data
M3 = left_join(mfr, q1, by="PREG_ID")
M3$FETID = paste(mfr$PREG_ID, mfr$CHILDNUM, sep="_")
M3 = M3[,c("FETID", "PREG_ID", "CHILDNUM", "AGE", "GA", "SEX", "BIRTHWEIGHT", "PARITY", "AA85", "AA86", "AA87")]
names(M3)[9:11] = c("WEIGHT", "WEIGHT_1stT", "HEIGHT")
head(M3); dim(M3); length(unique(M3$FETID))

# initial cleanup

## abnormal values should be removed
## (although abnormal SEX values usually are not present after filtering malformations etc.)
M3[(M3$SEX<1 | M3$SEX>2), "SEX"] = NA
M3$GA<-M3$GA/7
M3[M3$GA<15, "GA"] = NA
M3[M3$BIRTHWEIGHT==0,"BIRTHWEIGHT"] = NA

# check for missing or strange values
table(M3$PARITY,useNA="a")
table(round(M3$GA,0),useNA="a")
table(M3$SEX,useNA="a")
table(round(M3$BIRTHWEIGHT,-2),useNA="a")

#####################
### CASE 1: Skjaerven

## based on:
# Skjaerven et al. (2001) "Birthweight by gestational age in Norway"
# http://onlinelibrary.wiley.com/doi/10.1034/j.1600-0412.2000.079006440.x/abstract

## read in the reference data table
ref = read.table(skjaervenfile,header=T)
rboys = ref[ref$SEX==1,]
rgirls = ref[ref$SEX==2,]

## given the GA and SEX of each fetus, assign the mean weight from reference growth curves
mboys=(M3$SEX==1 & !is.na(M3$GA))
mgirls=(M3$SEX==2 & !is.na(M3$GA))
mmiss=(is.na(M3$SEX) & !is.na(M3$GA))
M3$REFMEAN=rep(0, nrow(M3))
M3$REFSD=rep(0, nrow(M3))

### Setting rule=2 in the approx functions would make it use extreme endpoints instead of generating NA,
### when the values fall outside the range. Not recommended.

## assign means
M3$REFMEAN[which(mboys)] = approx(rboys$GA, rboys$MEANWEIGHT, M3$GA[which(mboys)])$y
M3$REFMEAN[which(mgirls)] = approx(rgirls$GA, rgirls$MEANWEIGHT, M3$GA[which(mgirls)])$y
### when SEX is NA, just use both genders and take the average
M3$REFMEAN[which(mmiss)] = (approx(rgirls$GA, rgirls$MEANWEIGHT, M3$GA[mmiss])$y +
                                    approx(rboys$GA, rboys$MEANWEIGHT, M3$GA[which(mmiss)])$y)/2

## assign SDs
M3$REFSD[which(mboys)]=approx(rboys$GA, rboys$SD, M3$GA[which(mboys)])$y
M3$REFSD[which(mgirls)]=approx(rgirls$GA, rgirls$SD, M3$GA[which(mgirls)])$y
M3$REFSD[which(mmiss)]=(approx(rgirls$GA, rgirls$SD, M3$GA[which(mmiss)])$y +
                                approx(rboys$GA, rboys$SD, M3$GA[which(mmiss)])$y)/2

## for the output, calculate the percentile of weight, using reference MEAN and SD
M3$PCTskjaerven<-pnorm(M3$BIRTHWEIGHT,M3$REFMEAN,M3$REFSD)*100
## SGA can then be determined simply
M3$SGAskjaerven<-as.numeric(M3$PCTskjaerven<10)


#####################
### CASE 2: HADLOCK

## based on:
# Hadlock FP et al. (1991) "In utero analysis of fetal growth: a sonographic weight standard",
# with correction for GA < 25 weeks as described in:
# http://www.gestation.net/GROW_documentation.pdf, page 5

hadlock = function(ga,weight,refmean,refsd){
        ## input requires GA in rounded weeks
        ## weight = birthweight (in g)
        ## sex = 1 (boys) or 2 (girls)
        ## refmean and refsd are calculated at 280 days

        ## original Hadlock equation is valid only for GA >=25
        prop25=(299.1 - 31.85 * ga + 1.094 * ga^2 - 0.01055 * ga^3)/100
        prop24=(-5.86048381 + 1.419180433 * ga - 0.116517911 * ga^2 + 0.004154453 * ga^3)/100
        prop=ifelse(ga>=25, prop25, prop24)
        
        ## limits for reasonable prediction
        prop[which(ga>42.5 | ga<19.5)]=NA
        
        pnorm(weight, prop*refmean, prop*refsd)*100
}

## calculate parameters @ 280 days
mw40_boys = mean(M3[which(M3$GA==40 & M3$SEX==1),"BIRTHWEIGHT"], na.rm=T)
mw40_girls = mean(M3[which(M3$GA==40 & M3$SEX==2),"BIRTHWEIGHT"], na.rm=T)
sd40_boys = sd(M3[which(M3$GA==40 & M3$SEX==1),"BIRTHWEIGHT"], na.rm=T)
sd40_girls = sd(M3[which(M3$GA==40 & M3$SEX==2),"BIRTHWEIGHT"], na.rm=T)

## instead of Skjaerven means, assign the above constants to the df for more flexible calculations
M3$REFMEAN[which(mboys)]=mw40_boys
M3$REFMEAN[which(mgirls)]=mw40_girls
M3$REFMEAN[which(mmiss)]=mean(mw40_boys,mw40_girls)
M3$REFSD[which(mboys)]=sd40_boys
M3$REFSD[which(mgirls)]=sd40_girls
M3$REFSD[which(mmiss)]=mean(sd40_boys,sd40_girls)

## for the output, calculate the percentile of weight, using reference MEAN and SD
M3$PCThadlock<-hadlock(M3$GA,M3$BIRTHWEIGHT,M3$REFMEAN,M3$REFSD)
## SGA can then be determined simply
M3$SGAhadlock<-as.numeric(M3$PCThadlock<10)


#####################
### CASE 3: GARDOSI

## parameters taken from:
# Gardosi, Clausson & Francis (2009) "The value of customised centiles in
# assessing perinatal mortality risk associated with parity and maternal size"
# http://onlinelibrary.wiley.com/doi/10.1111/j.1471-0528.2009.02245.x/full

## for the original description of the method, see:
# Gardosi et al. (1995) "An adjustable fetal weight standard"
# http://onlinelibrary.wiley.com/doi/10.1046/j.1469-0705.1995.06030168.x/pdf

gardosi = function(weight,height,sex,parity,ga,birthweight) {
        ## sex       sex of the fetus (coded 2 for girls, 1 for boys)
        ## parity    parity (first child is 0)
        ## ga        gestational age in rounded weeks
        
        ## note that the means of mother height and weight do not match ours - neither does the sex correction -
        ## but the resulting predictions do not show any systemic bias
        
        w = weight - 65    # mean center (1st trim weight)
        wa = 9.066 * w - 0.067 * w^2
        h = height - 166      # mean center
        ha = 8.316*h - 0.006 * h^3
        paritycorr = rep(20,189)        # limited to 20 kids, but should be more than enough
        paritycorr[1:4] = c(0,136,174.4,183.4)
        p = paritycorr[parity+1]
        s = (-1)^(sex+1) * 64.1
        
        tow = 3575.2 + wa + ha + s + p          # predicted weight at 40 weeks GA
        tsd = tow * 0.11        # 11 % coef. of var. from original article (gardosi 1995) also matches our data
        hadlock(ga, birthweight, tow, tsd)
}

M3$PCTgardosi<-gardosi(M3$WEIGHT_1stT,M3$HEIGHT,M3$SEX,M3$PARITY,M3$GA,M3$BIRTHWEIGHT)
M3$SGAgardosi<-as.numeric(M3$PCTgardosi<10)

#####################
### CASE 4: MARSAL

## based on:
# Marsal et al. (1996) "Intrauterine growth curves based on ultrasonically estimated foetal weights"
# http://onlinelibrary.wiley.com/doi/10.1111/j.1651-2227.1996.tb14164.x/pdf

marsal = function(ga,birthweight,sex){
        if(is.na(sex)) sex=0
        if(sex==1){
                ## boys
                mw=-1.907345e-6*ga^4 + 1.140644e-3*ga^3 - 1.336265e-1*ga^2 + 1.976961*ga + 2.410053e+2
        } else if(sex==2){
                ## girls
                mw=-2.761948e-6*ga^4 + 1.744841e-3*ga^3 - 2.893626e-1*ga^2 + 18.91197*ga - 4.135122e+2
        } else {
                ## ?
                mw=-2.278843e-6*ga^4 + 1.402168e-3*ga^3 - 2.008726e-1*ga^2 + 9.284121*ga - 4.125956e+1
        }
        pnorm(birthweight, mw, abs(0.12*mw))*100
}

M3$PCTmarsal=sapply(1:nrow(M3), function(x) marsal(M3$GA[x]*7,M3$BIRTHWEIGHT[x],M3$SEX[x]))
## note that the cutoff here is -2 SDs, not 10th percentile
M3$SGAmarsal=as.numeric(M3$PCTmarsal<pnorm(-2)*100)

#####################
### CASE 5: EMPIRICAL PERCENTILES

## simply converts birthweight to a percentile score
## of all weights with the same GA and SEX.
#### !! note that NA birthweights are removed here, as they interfere with percentile calculations !!
M4 = group_by(M3, GA, SEX) %>% filter(!is.na(BIRTHWEIGHT))
M4 = mutate(M4, PCTempirical = dense_rank(BIRTHWEIGHT)) %>%
        mutate(PCTempirical = (PCTempirical-0.5)/max(PCTempirical)*100)
M4$SGAempirical<-as.numeric(M4$PCTempirical<10)


####################################################################################
####################################################################################

##### RESULTS #####

par(mfrow=c(1,2))

## check the distributions; red line indicates the SGA cutoff
hist(M4$PCTskjaerven, breaks=100, col="grey", main = "Method: Skjaerven")
abline(v=10, col="red")
text(paste("total number of SGAs:", sum(M4$SGAskjaerven,na.rm=T)), x=50, y=1.5*nrow(M4)/100, col="darkred")

hist(M4$PCThadlock, breaks=100, col="grey", main = "Method: Hadlock")
abline(v=10, col="red")
text(paste("total number of SGAs:", sum(M4$SGAhadlock,na.rm=T)), x=50, y=1.5*nrow(M4)/100, col="darkred")

hist(M4$PCTgardosi, breaks=100, col="grey", main = "Method: Gardosi")
abline(v=10, col="red")
text(paste("total number of SGAs:", sum(M4$SGAgardosi,na.rm=T)), x=50, y=1.5*nrow(M4)/100, col="darkred")

hist(M4$PCTmarsal, breaks=100, col="grey", main = "Method: Marsal")
abline(v=10, col="red")
text(paste("total number of SGAs:", sum(M4$SGAmarsal,na.rm=T)), x=50, y=1.5*nrow(M4)/100, col="darkred")

hist(M4$PCTempirical, breaks=100, col="grey", main = "Method: Empirical")
abline(v=10, col="red")
text(paste("total number of SGAs:", sum(M4$SGAempirical,na.rm=T)), x=50, y=1.5*nrow(M4)/100, col="darkred")

## inspect the differences in SGA indications
### note that tests can produce NAs if the GA value falls outside normal range
table(M3$SGAskjaerven, M3$SGAhadlock, dnn=c("Skjaerven","Hadlock"), useNA="a")
table(M3$SGAskjaerven, M3$SGAgardosi, dnn=c("Skjaerven","Gardosi"), useNA="a")
table(M3$SGAskjaerven, M3$SGAmarsal, dnn=c("Skjaerven","Marsal"), useNA="a")
table(M3$SGAhadlock, M3$SGAgardosi, dnn=c("Hadlock","Gardosi"), useNA="a")
table(M3$SGAhadlock, M3$SGAmarsal, dnn=c("Hadlock","Marsal"), useNA="a")
table(M3$SGAgardosi, M3$SGAmarsal, dnn=c("Gardosi","Marsal"), useNA="a")
table(M4$SGAskjaerven, M4$SGAempirical, dnn=c("Skjaerven","Empirical"), useNA="a")
table(M4$SGAhadlock, M4$SGAempirical, dnn=c("Skjaerven","Empirical"), useNA="a")
table(M4$SGAgardosi, M4$SGAempirical, dnn=c("Skjaerven","Empirical"), useNA="a")
table(M4$SGAmarsal, M4$SGAempirical, dnn=c("Skjaerven","Empirical"), useNA="a")

par(mfrow=c(1,1))
venn(na.exclude(M3[, grep("^SGA", colnames(M3))]))

## write final output (fetal info + birthweight percentiles)
out=M4[,c("PREG_ID","CHILDNUM","GA","SEX","BIRTHWEIGHT",
          "PCTskjaerven","PCThadlock","PCTgardosi","PCTmarsal","PCTempirical",
          "SGAskjaerven","SGAhadlock","SGAgardosi","SGAmarsal","SGAempirical")]
## converts GA back to days
out$GA=out$GA*7
head(out); dim(out)

write.table(format(out, digits=4), file_out, sep="\t", quote=F,row.names=F,col.names=T)
