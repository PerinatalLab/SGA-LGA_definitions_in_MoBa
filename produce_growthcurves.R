
#### This script estimates normal growth curves from population data on healthy births.
#### Input data, GOODBIRTHS set, has been generated from MFR, as described below:

## SQL query used to select good births:
# select PREG_ID_540,SVLEN_UL_DG,VEKT,LENGDE,HODE,KJONN from mfr where
# DODFODTE_5=0 AND SPABORT_12_5=0 AND SPABORT_23_5=0   #### this is omitted in goodbirths2!
# AND SVLEN_UL_DG>0 AND FSTART=1 AND KSNITT=0 AND FLERFODSEL=0
# AND (TILSTAND=5 OR TILSTAND<3) AND C00_MALF_ALL=0 AND IVF=0
# INTO OUTFILE '/tmp/output.csv' FIELDS TERMINATED BY '\t';

## meaning:
# no previous stillborns or miscarriages (omitted in goodbirths2!)
# only spontaneous birth (no induction/cesarean)
# only singletons
# only those born alive (including "died later" or "no follow-up")
# no malformations
# no ivf
# GA determined by ultrasound


library(dplyr)
library(ggplot2)

## read in skjaerven's table
ref=read.table("NSGA_BirthWeight_by_Skjaerven_2000_MEANSD.txt",header=T)
rboys<-ref[ref$SEX==1,]
rgirls<-ref[ref$SEX==2,]

## read in our data
dat=read.table("output_mfr_goodbirths2.csv",h=T)
head(dat); dim(dat)
dat=dat[which(dat$WEIGHT>0),]
## use TRUNC for Skjaerven (completed weeks), otherwise ROUND
#dat$GAV=trunc(dat$GA/7)
dat$GAV=round(dat$GA/7,0)

## remove nonsense
dat$WEIGHT[which(dat$WEIGHT==0)]=NA
dat$LENGTH[which(dat$LENGTH==0)]=NA
dat$HEADSIZE[which(dat$HEADSIZE==0)]=NA

## calculate per-week stats
dd=ddply(dat, c("GAV","SEX"), function(df) data.frame(m=mean(df$WEIGHT,na.rm=T), sd=sd(df$WEIGHT,na.rm=T),
                                q10=quantile(df$WEIGHT, 0.1, na.rm=T), q90=quantile(df$WEIGHT, 0.9, na.rm=T), n=nrow(df) ) )
dd2=dd[which(dd$n>1),]
ddboys=dd2[which(dd2$SEX==1),]
ddboys=merge(ddboys,rboys,by.x="GAV",by.y="GA",all.x=T)[,-c(2,8)]
names(ddboys)[7]="REFM"; names(ddboys)[8]="REFSD"
ddboys$REFQ10=qnorm(0.1,ddboys$REFM,ddboys$REFSD)

ddgirls=dd2[which(dd2$SEX==2),]
ddgirls=merge(ddgirls,rgirls,by.x="GAV",by.y="GA",all.x=T)[,-c(2,8)]
names(ddgirls)[7]="REFM"; names(ddgirls)[8]="REFSD"
ddgirls$REFQ10=qnorm(0.1,ddgirls$REFM,ddgirls$REFSD)

## gardosi things
mw40_boys=mean(dat[which(dat$GA==280 & dat$SEX==1),"WEIGHT"], na.rm=T)
mw40_girls=mean(dat[which(dat$GA==280 & dat$SEX==2),"WEIGHT"], na.rm=T)
sd40_boys=sd(dat[which(dat$GA==280 & dat$SEX==1),"WEIGHT"], na.rm=T)
sd40_girls=sd(dat[which(dat$GA==280 & dat$SEX==2),"WEIGHT"], na.rm=T)

gardosi=function(ga,meanweight){
        ## correction if this method uses true weeks, not completed ones
        #ga=ga+0.5
        prop=(298.8 - 31.85 * ga + 1.094 * ga^2 - 0.01055 * ga^3)/100
        prop*meanweight
}

## calculate the actual reference curves
ddboys$gardosi=gardosi(ddboys$GAV,mw40_boys)
ddgirls$gardosi=gardosi(ddgirls$GAV,mw40_girls)
ddboys$gardosi10=gardosi(ddboys$GAV,qnorm(0.1,mw40_boys,sd40_boys))
ddgirls$gardosi10=gardosi(ddgirls$GAV,qnorm(0.1,mw40_girls,sd40_girls))
ddboys$gardosi90=gardosi(ddboys$GAV,qnorm(0.9,mw40_boys,sd40_boys))
ddgirls$gardosi90=gardosi(ddgirls$GAV,qnorm(0.9,mw40_girls,sd40_girls))

## bind everything for plotting
gg=rbind(data.frame(ddboys,SEX="M"),data.frame(ddgirls,SEX="F"))
gg=gg[which(gg$n>2),]

## create the smoothed percentile lines for ribbon plot
p=ggplot(data=gg,aes(x=GAV)) + facet_wrap(~SEX) +
        geom_smooth(aes(y=q10),se=FALSE) + geom_smooth(aes(y=q90),se=FALSE)
plow=ggplot_build(p)$data[[1]]
pupp=ggplot_build(p)$data[[2]]
names(plow)=names(pupp)=c("GAV","min","SEX","group")
plow$max=pupp$min
plow$SEX=as.numeric(plow$SEX)
plow$SEX[which(plow$SEX==1)]="M"
plow$SEX[which(plow$SEX==2)]="F"

## using 1.282SDs
ggplot(data=gg,aes(x=GAV))+facet_wrap(~SEX)+labs(x="Gestational age, weeks",y="Birthweight, grams")+
        geom_point(aes(y=m),col="blue4")+geom_smooth(aes(y=m),se=FALSE,col="darkblue")+
        geom_line(aes(y=REFM),col="black",linetype=2)+
        geom_linerange(aes(ymin=qnorm(0.1,m,sd),ymax=qnorm(0.9,m,sd)), col="royalblue1")+
        geom_ribbon(data=plow,aes(x=GAV,ymin=min,ymax=max),fill="lightblue2",alpha=0.3)+
        geom_line(aes(y=qnorm(0.1,REFM,REFSD)),linetype=2,col="grey40")+
        geom_line(aes(y=qnorm(0.9,REFM,REFSD)),linetype=2,col="grey40")

## using actual 10-90 percentiles
ggplot(data=gg,aes(x=GAV))+facet_wrap(~SEX)+labs(x="Gestational age, weeks",y="Birthweight, grams")+
        geom_point(aes(y=m),col="blue4")+geom_smooth(aes(y=m),se=FALSE,col="darkblue")+
        geom_line(aes(y=REFM),col="black",linetype=2)+
        geom_linerange(aes(ymin=q10,ymax=q90), col="royalblue1")+
        geom_ribbon(data=plow,aes(x=GAV,ymin=min,ymax=max),fill="lightblue2",alpha=0.3)+
        geom_line(aes(y=qnorm(0.1,REFM,REFSD)),linetype=2,col="grey40")+
        geom_line(aes(y=qnorm(0.9,REFM,REFSD)),linetype=2,col="grey40")+
        geom_bar(stat="identity",aes(y=n/4.3),fill="grey60",alpha=0.4)+geom_text(aes(label=n,y=n/4.3+100),size=3)

## using actual 10-90 percentiles and gardosi prediction
ggplot(data=gg,aes(x=GAV))+facet_wrap(~SEX)+labs(x="Gestational age, weeks",y="Birthweight, grams")+
        geom_point(aes(y=m),col="blue4")+geom_smooth(aes(y=m),se=FALSE,col="darkblue")+
        geom_line(aes(y=gardosi),col="black",linetype=2)+
        geom_linerange(aes(ymin=q10,ymax=q90), col="royalblue1")+
        geom_ribbon(data=plow,aes(x=GAV,ymin=min,ymax=max),fill="lightblue2",alpha=0.3)+
        geom_line(aes(y=gardosi10),linetype=2,col="grey40")+
        geom_line(aes(y=gardosi90),linetype=2,col="grey40")+
        geom_bar(stat="identity",aes(y=n/4.3),fill="grey60",alpha=0.4)+geom_text(aes(label=n,y=n/4.3+100),size=3)

#####################
### we can do our own gardosi-thingy
inv_hadlock=function(ga,weight){
        # calculates the tow, given weight at any GA
        ga=ga/7
        prop=(298.8 - 31.85 * ga + 1.094 * ga^2 - 0.01055 * ga^3)/100
        weight/prop
}

normies=M2[which(M2$GA>258 & M2$GA<295),]
attach(normies)
rm(BIRTHWEIGHT, SEX, WEIGHT_1stT, HEIGHT, PARITY, W2, W3, H2, H3, P2, P3, P4)
BIRTHWEIGHT=inv_hadlock(GA, BIRTHWEIGHT)
SEX=SEX-1
WEIGHT_1stT=WEIGHT_1stT-mean(WEIGHT_1stT,na.rm=T)
W2=WEIGHT_1stT^2
W3=WEIGHT_1stT^3
HEIGHT=HEIGHT-mean(HEIGHT,na.rm=T)
H2=HEIGHT^2
H3=HEIGHT^3
P2=as.numeric(PARITY==2)
P3=as.numeric(PARITY==3)
P4=as.numeric(PARITY==4)
PARITY=as.numeric(PARITY==1)
nr=lm(BIRTHWEIGHT~SEX+WEIGHT_1stT+W2+W3+HEIGHT+H2+H3+PARITY+P2+P3+P4)
summary(nr)
detach(normies)

data.frame(our=coef(nr),theirs=c(3478.4,58.4*2,8.7,-)),