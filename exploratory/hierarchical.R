## --- libraries
library(here)
library(rstan)
library(data.table)
library(ggplot2)
library(readxl)
library(posterior)

## --- plotting function
source(here('exploratory/corplot.R')) #need lattice and hexbin installed

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
smps[,c('.chain','.iteration','.draw')] <- NULL

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
