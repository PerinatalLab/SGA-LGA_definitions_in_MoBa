####################################################################################
##### CALCULATE BIRTH-WEIGHT PERCENTILES (SKJAERVEN, HADLOCK, GARDOSI, MARSAL) #####

## This script is intended to be used and distributed together with
## Python GUI launcher defineSGA.exe.

### read input parameters

args = commandArgs(TRUE)
infile = args[1]
mfrfile = args[2]
file_dir = paste(args[3], "/")
outFileStem = args[4]
## assuming github is not used
## assume data tables for Skjaerven's method were in the same folder
skjaervenfile = "NSGA_BirthWeight_by_Skjaerven_2000_MEANSD.txt"

calculateSkj = (args[5]==1)
calculateHad = (args[6]==1)
calculateGar = (args[7]==1)
calculateMar = (args[8]==1)

includeMfrCols = args[9:13]
includeQ1Cols = args[14:16]
includePerc = (args[17]==1)
includeSGA = (args[18]==1)

infileSep = args[19]
mfrfileSep = args[20]

infileHead = (args[21]==1)
mfrfileHead = (args[22]==1)

##### INPUT REQUIREMENTS #####

## This script requires three input files:
##
## 1. imputed maternal data, in the following format:
##      (INCLUDING HEADER LINE) sep \t
##      PREG_ID     AA85    AA86    AA87    flAA85  flAA86  flAA87
##
## 2. additional information directly from MFR, in the following format:
##      (NO HEADER LINE) sep ,
##      PREG_ID     CHILD_NUMBER    MOTHERS_AGE     GA(days)      SEX(1/2)     BIRTHWEIGHT(g)     PARITY
##
## 3. population means and standard deviations for Skjaerven's method, in the following format:
##      (INCLUDING HEADER LINE) sep whitespace
##      GA(completed weeks)       SEX(1/2)       MEAN      SD

####################################################################################

##### CALCULATIONS #####
## install the SQLDF and GPLOTS packages if you do not have it yet
if(!require(sqldf)) install.packages('sqldf')
library(sqldf)
if(!require(gplots)){ install.packages('gplots') }
library(gplots)

## get the date stamp
date_stamp = paste(unlist(strsplit(substr(Sys.time(),1,10),"-")), collapse="")
file_out = paste(file_dir, outFileStem, "_", date_stamp, ".txt",sep="")
file_out_pdf = paste(file_dir, outFileStem, "_", date_stamp, ".pdf",sep="")

## read input files
q1 = read.csv(infile, sep=infileSep, header=infileHead)
names(q1)[1:4] = c("PREG_ID","AA85","AA86","AA87")
head(q1); dim(q1)

mfr = read.csv(mfrfile, sep=mfrfileSep, header=mfrfileHead)
names(mfr)[1:7] = c("PREG_ID","CHILDNUM","AGE","GA","SEX","BIRTHWEIGHT","PARITY")
head(mfr); dim(mfr)

## merge them
### note that a unique fetal ID is also created - it can be useful to keep track of the data
M3=sqldf("SELECT mfr.PREG_ID || '_' || mfr.CHILDNUM as FETID,
        mfr.PREG_ID, mfr.CHILDNUM, mfr.AGE, mfr.GA, mfr.SEX, mfr.BIRTHWEIGHT, mfr.PARITY,
        q1.AA85 as WEIGHT, q1.AA86 as WEIGHT_1stT, q1.AA87 as HEIGHT
        FROM mfr LEFT JOIN q1 ON mfr.PREG_ID=q1.PREG_ID")
head(M3); dim(M3); length(unique(M3$FETID))
M3=merge(mfr, q1, by="PREG_ID", all.x=TRUE)
head(M3)

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

if(calculateSkj){
        ## read in the reference data table
        ref = read.table(skjaervenfile, header=T)
        rboys = ref[ref$SEX==1,]
        rgirls = ref[ref$SEX==2,]
        
        ## given the GA and SEX of each fetus, assign the mean weight from reference growth curves
        mboys=(M3$SEX==1 & !is.na(M3$GA))
        mgirls=(M3$SEX==2 & !is.na(M3$GA))
        mmiss=(is.na(M3$SEX) & !is.na(M3$GA))
        M3$REFMEAN=NULL
        M3$REFSD=NULL
        
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
}

#####################
### CASE 2: HADLOCK

## based on:
# Hadlock FP et al. (1991) "In utero analysis of fetal growth: a sonographic weight standard",
# with correction for GA < 25 weeks as described in:
# http://www.gestation.net/GROW_documentation.pdf, page 5

if(calculateHad){
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
}

#####################
### CASE 3: GARDOSI

## parameters taken from:
# Gardosi, Clausson & Francis (2009) "The value of customised centiles in
# assessing perinatal mortality risk associated with parity and maternal size"
# http://onlinelibrary.wiley.com/doi/10.1111/j.1471-0528.2009.02245.x/full

## for the original description of the method, see:
# Gardosi et al. (1995) "An adjustable fetal weight standard"
# http://onlinelibrary.wiley.com/doi/10.1046/j.1469-0705.1995.06030168.x/pdf

if(calculateGar){
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
}

#####################
### CASE 4: MARSAL

## based on:
# Marsal et al. (1996) "Intrauterine growth curves based on ultrasonically estimated foetal weights"
# http://onlinelibrary.wiley.com/doi/10.1111/j.1651-2227.1996.tb14164.x/pdf

if(calculateMar){
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
}

####################################################################################

##### RESULTS #####

## restore GA days
M3$GA = M3$GA*7
## default columns
out = M3[,c("PREG_ID","CHILDNUM")]
## add age, ga, sex, birthweight, parity
out = cbind(out, M3[,which(includeMfrCols==1)+3])
## add weight, weight_1stT, height
out = cbind(out, M3[,which(includeQ1Cols==1)+8])

## add percentile/SGA indicator columns and plot outputs into a pdf
pdf(file_out_pdf)
if(calculateSkj){
        if(includePerc){ out = cbind(out, M3$PCTskjaerven) }
        if(includeSGA){ out = cbind(out, M3$SGAskjaerven) }
        
        hist(M3$PCTskjaerven, breaks=100, col="grey", main = "Method: Skjaerven", xlab="Weight percentile")
        abline(v=10, col="red")
        text(paste("total number of SGAs:", sum(M3$SGAskjaerven,na.rm=T), "/", nrow(M3)), x=50, y=1.5*nrow(M3)/100, col="darkred")
}
if(calculateHad){
        if(includePerc){ out = cbind(out, M3$PCThadlock) }
        if(includeSGA){ out = cbind(out, M3$SGAhadlock) }
        
        hist(M3$PCThadlock, breaks=100, col="grey", main = "Method: Hadlock", xlab="Weight percentile")
        abline(v=10, col="red")
        text(paste("total number of SGAs:", sum(M3$SGAhadlock,na.rm=T), "/", nrow(M3)), x=50, y=1.5*nrow(M3)/100, col="darkred")
}
if(calculateGar){
        if(includePerc){ out = cbind(out, M3$PCTgardosi) }
        if(includeSGA){ out = cbind(out, M3$SGAgardosi) }
        
        hist(M3$PCTgardosi, breaks=100, col="grey", main = "Method: Gardosi", xlab="Weight percentile")
        abline(v=10, col="red")
        text(paste("total number of SGAs:", sum(M3$SGAgardosi,na.rm=T), "/", nrow(M3)), x=50, y=1.5*nrow(M3)/100, col="darkred")
}
if(calculateMar){
        if(includePerc){ out = cbind(out, M3$PCTmarsal) }
        if(includeSGA){ out = cbind(out, M3$SGAmarsal) }
        
        hist(M3$PCTmarsal, breaks=100, col="grey", main = "Method: Marsal", xlab="Weight percentile")
        abline(v=10, col="red")
        text(paste("total number of SGAs:", sum(M3$SGAmarsal,na.rm=T), "/", nrow(M3)), x=50, y=1.5*nrow(M3)/100, col="darkred")
}

## inspect the differences in SGA indications
venn(M3[, grep("^SGA", colnames(M3))])
dev.off()

write.table(format(out, digits=4), file_out, sep="\t", quote=F,row.names=F,col.names=T)
