## --- libraries
library(here)
library(rstan)
library(data.table)
library(ggplot2)
library(readxl)
library(posterior)
library(corplot) #devtools::install_github('petedodd/corplot')
library(stringr)
library(ggrepel)
library(patchwork)

## --- compile stan model
fn <- here('exploratory/stan/hier.stan')
hm <- stan_model(file=fn)

## --- read data
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


## --- make data for stan
tmp <- dcast(W,ward ~ year,value.var='tb')
tmp[,ward:=NULL]
tmp <- as.matrix(tmp)
str(tmp)

## list
sdata <- list(
  K = nrow(tmp),
  NT = ncol(tmp),
  tACF = as.integer(which(colnames(tmp)=='1957')),
  Y = tmp
)

## --- sample from stan model
n <- 5
samps <- sampling(hm,data=sdata,chains=n,cores=n,iter=2e3)

## --- look at posteriors
smps <- as_draws_df(samps)
smps <- mutate_variables(smps, RR.peak = exp(`gm[1]`))
smps <- mutate_variables(smps, RR.level = exp(`gm[2]`))
smps <- mutate_variables(smps, RR.slope = exp(`gm[3]`))
smps <- subset_draws(smps,variable=c('RR.peak','RR.level','RR.slope'))
smpsf <- smps #with the chain/iteration/draw included
drp <- c('.chain','.iteration','.draw')
smps[,drp] <- NULL

fn <- here('exploratory/figures/post.real.pdf')
corplot(smps,file=fn)


emps <- subset_draws(as_draws_df(samps),variable='gm')
emps[,c('.chain','.iteration','.draw')] <- NULL

fn <- here('exploratory/figures/post.log.pdf')
corplot(emps,file=fn)

## --- posterior summaries
(dtmp <- summarise_draws(smpsf))
fn <- here('exploratory/figures/post.real.summary.csv')
fwrite(dtmp,file=fn)


## --- site level effects
smps <- subset_draws(smps,variable=c('RR.peak','RR.level','RR.slope'))

smps <- as_draws_df(samps)
smpsr <- subset_draws(smps,variable='RRs')
smpsrm <- melt(as.data.table(smpsr),id=drp)

## functions for reformatting outputs
getsno <- function(x){
  a <- str_extract(x,"\\[(.*?),")
  a <- gsub('\\[','',a)
  a <- gsub(',','',a)
  as.integer(a)
}
getrno <- function(x){
  a <- str_extract(x,",(.*?)\\]")
  a <- gsub('\\]','',a)
  a <- gsub(',','',a)
  as.integer(a)
}
## keys
rcts <- c('RR.peak','RR.level','RR.slope')
tmp2 <- dcast(W,ward ~ year,value.var='tb')
scts <- as.character(tmp2$ward)
getw <- function(x) scts[getsno(x)]
getv <- function(x) rcts[getrno(x)]
## test
getw(head(smpsrm$variable))
getv(head(smpsrm$variable))

## apply to data
smpsrm[,ward:=getw(variable)]
smpsrm[,effect:=getv(variable)]

peakvlvl <- dcast(smpsrm[effect != 'RR.slope'],.draw + ward ~ effect,value.var = 'value')
peakvlvlm <- peakvlvl[,.(RR.level=mean(RR.level),RR.peak=mean(RR.peak)),by=ward]

## NOTE each point +ve correlation
ggplot(peakvlvl,aes(RR.level,RR.peak,group=ward))+
  stat_ellipse()+
  geom_point(data=peakvlvlm,col=2)

## pairs for means:
smy <- smpsrm[,.(value=mean(value)),by=.(ward,effect)]
smy <- dcast(smy,ward  ~ effect,value.var='value')

## --- plot of ward-level mean effects
gg1 <- ggplot(smy,aes(RR.level,RR.peak,label=ward))+
  geom_point()+
  geom_text_repel()+
  ggtitle('Level vs. Peak')
gg2 <- ggplot(smy,aes(RR.level,RR.slope,label=ward))+
  geom_point()+
  geom_text_repel()+
  ggtitle('Level vs. Slope')
gg3 <- ggplot(smy,aes(RR.peak,RR.slope,label=ward))+
  geom_point()+
  geom_text_repel()+
  ggtitle('Peak vs. Slope')


gg <- gg1+gg2+gg3
ggsave(gg,file=here('exploratory/figures/post.ward.meanRRs.pdf'),h=5,w=15)
