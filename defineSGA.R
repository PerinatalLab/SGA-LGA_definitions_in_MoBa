########################################################################
####### CONVERT BIRTH-WEIGHT TO PERCENTILES  (SKJAERVEN, HADLOCK, GARDOSI, MARSAL)

#### This script uses four different methods to evaluate birth-weight.
#### The output is a percentile and SGA indication (1/0) for each fetus.

# read in the csv data files
## (note the required columns and their order)

## get the date stamp
date_stamp = paste(unlist(strsplit(substr(Sys.time(),1,10),"-")),collapse="")
## get the hash of the current state of the Git folder
setwd("~/Documents/gitrep/SGA-LGA_definitions_in_MoBa/")
hash = system(paste("git log --pretty=format:'%h' -n 1"),intern=TRUE)

file_dir = "~/Desktop/MoBa_v6/"
file_out = paste(file_dir,"MOBA_PDB1581_fetalWEIGHTandSGA_",date_stamp,"_",hash,".txt",sep="")

## this file contains the cleaned, corrected and imputed mother height/weight info
#### (don't forget to set the INHASH each time a new data file is produced)
infile="MOBA_PDB1581_IMPUTED_maternalHgh1Wgh1Wgh2_20151005_"
INHASH="1a447dc"
q1=read.csv(paste(file_dir, infile, INHASH, ".txt", sep=""), sep="\t", header=T)
names(q1)=c("PREG_ID","AA85","AA86","AA87","flAA85","flAA86","flAA87")
head(q1); dim(q1)

## these two files contain info retrieved directly from MoBa without any changes
en=read.csv(paste(file_dir,"output_q2_kilojoules.csv",sep=""), sep=",", header=T)
names(en)=c("PREG_ID","KJ")
head(en); dim(en)

mfr=read.csv(paste(file_dir,"output_mfr_basicfetalinfo.csv",sep=""),sep=",",header=F)
names(mfr)=c("PREG_ID","CHILDNUM","AGE","GA","SEX","BIRTHWEIGHT","PARITY")
head(mfr); dim(mfr)

## merge the last two csvs
library(sqldf)
M3=sqldf("SELECT mfr.PREG_ID || '_' || mfr.CHILDNUM as FETID, mfr.*, q1.AA85 as WEIGHT, q1.AA86 as WEIGHT_1stT, q1.AA87 as HEIGHT,
        q1.flAA85, q1.flAA86, q1.flAA87 FROM mfr INNER JOIN q1 ON mfr.PREG_ID=q1.PREG_ID")
head(M3); dim(M3); length(unique(M3$FETID))

# initial cleanup

## abnormal values should be removed
## (although abnormal SEX values are not present in GOODBIRTH set)
M3[(M3$SEX<1 | M3$SEX>2), "SEX"] = NA
M3$GA<-M3$GA/7
M3[M3$GA<15, "GA"] = NA
M3[M3$BIRTHWEIGHT==0,"BIRTHWEIGHT"] = NA

# check for missing or strange values
table(M3$PARITY,useNA="a")
table(round(M3$GA,0),useNA="a")
table(M3$SEX,useNA="a")
table(round(M3$BIRTHWEIGHT,-2),useNA="a")

### CASE 1: skjaerven

## read in the table
ref=read.table("NSGA_BirthWeight_by_Skjaerven_2000_MEANSD.txt",header=T)
head(ref)
rboys<-ref[ref$SEX==1,]
rgirls<-ref[ref$SEX==2,]

## given the GA and SEX of each fetus, retrieve the mean weight from the reference growth curves
mboys=(M3$SEX==1 & !is.na(M3$GA))
mgirls=(M3$SEX==2 & !is.na(M3$GA))
mmiss=(is.na(M3$SEX) & !is.na(M3$GA))

### setting rule=2 would tell approx to use extreme endpoints instead of generating NA, when the values fall outside the range
### (don't think that's a good idea)
M3$REFMEAN=NULL
M3$REFSD=NULL
M3$REFMEAN[which(mboys)] = approx(rboys$GA, rboys$MEANWEIGHT, M3$GA[which(mboys)])$y
M3$REFMEAN[which(mgirls)] = approx(rgirls$GA, rgirls$MEANWEIGHT, M3$GA[which(mgirls)])$y
### when SEX is NA, but GA is known, just use both genders and take the average
M3$REFMEAN[which(mmiss)] = (approx(rgirls$GA, rgirls$MEANWEIGHT, M3$GA[mmiss])$y +
                                    approx(rboys$GA, rboys$MEANWEIGHT, M3$GA[which(mmiss)])$y)/2

## retrieve SD in the same way
M3$REFSD[which(mboys)]=approx(rboys$GA, rboys$SD, M3$GA[which(mboys)])$y
M3$REFSD[which(mgirls)]=approx(rgirls$GA, rgirls$SD, M3$GA[which(mgirls)])$y
M3$REFSD[which(mmiss)]=(approx(rgirls$GA, rgirls$SD, M3$GA[which(mmiss)])$y +
                                approx(rboys$GA, rboys$SD, M3$GA[which(mmiss)])$y)/2

## for the output, calculate the percentile of weight, using reference MEAN and SD
M3$PCTskjaerven<-pnorm(M3$BIRTHWEIGHT,M3$REFMEAN,M3$REFSD)*100
## SGA can then be determined simply
M3$SGAskjaerven<-as.numeric(M3$PCTskjaerven<10)

## check the distribution
hist(M3$PCTskjaerven,breaks=100,col="grey")
abline(v=10,col="red")
text(paste("total number of SGAs:", sum(M3$SGAskjaerven,na.rm=T)),x=40, y=1300,col="darkred")

### CASE 2: HADLOCK

hadlock=function(ga,weight,refmean,refsd){
        ## input requires GA in rounded weeks
        ## weight = birthweight (in g)
        ## sex = 1 (boys) or 2 (girls)
        ## refmean and refsd are calculated at 280 days
        
        ## use this correction if input is in completed weeks:
        # ga=ga+0.5
        
        ## original Hadlock equation is valid only for GA >=25
        prop25=(299.1 - 31.85 * ga + 1.094 * ga^2 - 0.01055 * ga^3)/100
        prop24=(-5.86048381 + 1.419180433 * ga - 0.116517911 * ga^2 + 0.004154453 * ga^3)/100
        prop=ifelse(ga>=25, prop25, prop24)
        
        ## limits for reasonable prediction
        prop[which(ga>42.5 | ga<19.5)]=NA
        
        pnorm(weight, prop*refmean, prop*refsd)*100
}

## calculate parameters @ 280 days
mw40_boys=mean(M3[which(M3$GA==40 & M3$SEX==1),"BIRTHWEIGHT"], na.rm=T)
mw40_girls=mean(M3[which(M3$GA==40 & M3$SEX==2),"BIRTHWEIGHT"], na.rm=T)
sd40_boys=sd(M3[which(M3$GA==40 & M3$SEX==1),"BIRTHWEIGHT"], na.rm=T)
sd40_girls=sd(M3[which(M3$GA==40 & M3$SEX==2),"BIRTHWEIGHT"], na.rm=T)

## instead of Skjaerven means, assign the constants to the df for more flexible calculations
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

## check if both tests indicate the same SGAs
## note that Skjaerven produces NAs when the GA value falls outside normal range
table(M3$SGAskjaerven,M3$SGAhadlock,dnn=c("Skjaerven","Hadlock"),useNA="a")


### CASE 3: GARDOSI

gardosi<-function(weight,height,sex,parity,ga,birthweight) {
        # sex       sex of the fetus (coded 2 for girls, 1 for boys)
        # parity    parity              first child is 0
        # ga        gestational age in weeks
        
        ## This algorithm uses mother information to produce a personalized expected weight at term.
        ## The coefficients used here come from article Gardosi 2009, "The value of customised centiles..."
        ## note that the means of mother height and weight do not match ours - neither does the sex correction
        ## but the resulting predictions do not show any systemic bias
        
        w = weight - 65    # mean center (1st trim)          
        wa = 9.066 * w - 0.067 * w^2
        h = height - 166      # mean center
        ha = 8.316*h - 0.006 * h^3
        paritycorr = rep(20,189)                 # Guess no one has more than 20 kids!
        paritycorr[1:4] = c(0,136,174.4,183.4)
        p = paritycorr[parity+1]
        s = (-1)^(sex+1) * 64.1
        
        tow = 3575.2 + wa + ha + s + p          # predicted weight at 40 weeks GA
        tsd = tow*0.11        # 11 % coef. of var. from original article (gardosi 1995) matches our data
        hadlock(ga, birthweight, tow, tsd)
}

M3$PCTgardosi<-gardosi(M3$WEIGHT_1stT,M3$HEIGHT,M3$SEX,M3$PARITY,M3$GA,M3$BIRTHWEIGHT)
M3$SGAgardosi<-as.numeric(M3$PCTgardosi<10)

### CASE 4: MARSAL

marsal=function(ga,birthweight,sex){
        ## equations from Marsal 1996, "intrauterine growth curves..."
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
M3$SGAmarsal=as.numeric(M3$PCTmarsal<pnorm(-2))

## write final output (per-fetus info + birthweight percentiles)
out=M3[,c("PREG_ID","CHILDNUM","GA","SEX","BIRTHWEIGHT",
          "PCTskjaerven","PCThadlock","PCTgardosi","PCTmarsal", "SGAskjaerven","SGAhadlock","SGAgardosi","SGAmarsal")]
out$GA=out$GA*7
head(out); dim(out)

write.table(format(out, digits=4), file_out, sep="\t", quote=F,row.names=F,col.names=T)
