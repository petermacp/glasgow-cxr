## libraries
library(data.table)
library(ggplot2)
library(readxl)
library(here)
library(GGally)
library(ggrepel)

## ===== exploratory data analysis: plots and empirical patterns

## popdensity etc
warcov <- readRDS(here('populations/ward_covariates.rds'))
warcov <- as.data.table(unique(warcov[,c('division','ward','area','total_population','year')]))
warcov <- warcov[year==1950] #use baseline
warcov <- warcov[,.(division,ward,pdens=1e3*total_population/area)]

## filename
xfn <- here("2023-11-28_glasgow-acf.xlsx")

## ward populations
W_pops <- as.data.table(read_xlsx(path = xfn, sheet = "ward_population"))
WP <- as.data.table(W_pops[,.(ward,year,total_population)])

##list all the sheets
all_sheets <- excel_sheets(xfn)

##get the ward tb sheets
shts <- grep('by_ward',all_sheets,value=TRUE)
L <- list()
for(nm in shts) L[[nm]] <- as.data.table(read_xlsx(xfn,sheet=nm))
L <- rbindlist(L)

## ward level counts
W <- L[,.(tb=sum(cases,na.rm=TRUE)),by=.(year,ward)]
W <- merge(W,WP,by=c('ward','year'))

## look: 1957 = intervention
## numbers
ggplot(W,aes(year,tb))+
  geom_point()+
  facet_wrap(~ward,ncol=6)+
  theme_linedraw()

## rates
ggplot(W,aes(year,1e5*tb/total_population))+
  geom_point()+
  facet_wrap(~ward,ncol=6)+
  theme_linedraw()

ggsave(file=here('exploratory/figures/data.png'),w=12,h=9)

## code periods
W[,period:=ifelse(year<1957,'before','after')]
W[year==1957,period:='during']

## mean rates
WM <- W[,.(rate=mean(1e5*tb/total_population)),by=.(ward,period)]
WM <- dcast(WM,ward ~ period,value.var='rate')
WM[,reldiff:=(before-after)/before] # notification drop relative to 'before'
WM[,relpeak:=during/before]         # peak notifications relative to 'before'
WM[,RR.level:=(after)/before] # notification drop relative to 'before'
WM[,absdiff:=before-after]

## 76% cov:
WM[,prev:=before * (relpeak-1)/0.76] #assume 76% coverage
WM[,CDR:=before * 3 / prev] #P = IT
WM[,CDR:=(CDR)/(1+CDR)] #CDR estimate (inverse odds)
WM[,Ibefore:=before/CDR]

## sensitivity 70% cov:
WM[,prev70:=before * (relpeak-1)/0.7] #assume 70% coverage
WM[,CDR70:=before * 3 / prev70] #P = IT
WM[,CDR70:=(CDR70)/(1+CDR70)] #CDR estimate (inverse odds)
WM[,Ibefore70:=before/CDR70]


## merge in covs
WM <- merge(WM,warcov,by='ward')

## larger correlation
## removed: relpeak out as ~ CDR, RR.level ~ reldiff
GP <- ggpairs(WM[,.(CDR,RR.level,Ipre=Ibefore,Npre=before,prev,pdens)]) 
GP +
  ggdist::theme_ggdist() +
  theme(panel.background = element_rect(fill=NA, colour="grey78"))

ggsave(GP,file=here('exploratory/figures/CDR.pairs.png'),w=8,h=8)


## for text:
WM[,summary(CDR)]


## larger correlation
## removed: relpeak out as ~ CDR, RR.level ~ reldiff
GP <- ggpairs(WM[,.(CDR=CDR70,RR.level,Ipre=Ibefore70,Npre=before,prev=prev70,pdens)]) 
GP

ggsave(GP,file=here('exploratory/figures/CDR.pairs70.png'),w=8,h=8)

WM[,summary(CDR70)]
