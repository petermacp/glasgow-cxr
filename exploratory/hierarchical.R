## libraries
library(here)
library(rstan)

fn <- here('exploratory/stan/hier.stan')
hm <- stan_model(file=fn)

## TODO below here is just for filenames
library(data.table)
library(ggplot2)
library(readxl)
library(here)
library(GGally)

## ===== exploratory data analysis: plots and empirical patterns

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

## quite confusing
ggpairs(WM[,.(before,reldiff,relpeak)])

ggsave(file=here('exploratory/figures/effect.pairs.png'),w=5,h=5)


## controlling for each
mdl <- lm(data=WM[,.(before,reldiff,relpeak)],reldiff ~ before + relpeak)
summary(mdl) #bigger relative diff for smaller relative peak (controlling for before)
## NOTE perhaps the opposite from what one would expect if more success implies bigger impact
## should probably control for trend as this could generate this?
## or, does it represent confounding with better off or better access areas?
## TODO

## ===== simple calculation of notifications averted with delta method
WT <- W[,.(tb=sum(tb,na.rm=TRUE),Pop=sum(total_population)),by=.(period,year)]
## WT[period!='during',t:=1:(nrow(WT)-1)] #time with 'during' omitted
WT[,t:=1:nrow(WT)] #time as integer

ggplot(WT,aes(year,tb))+geom_point()+theme_linedraw()
ggplot(WT,aes(t,tb))+geom_point()+theme_linedraw()


## Model Without control
WT[,Iafter:=ifelse(period=='after',1,0)] #switch

## linear model
mod.woc <- glm(data=WT[period!='during'],
               family=poisson,
               tb ~
                 offset(log(Pop))+
                 (1 + t) +           # k_1 + s_1.t
                 Iafter + Iafter:t   # I(t)(a+b.t)
               )

summary(mod.woc)


## Define the functions for delta method (see Rachael Hit TB work)
D.fun <- function(theta){
  list2env(theta,envir = environment())
  vc <- Pop*(exp(a+k+(s+b)*itz)-exp(k+s*itz))
  -sum(vc[itz>intT]) #restrict to post intervention
}

## Need to differentiate by hand as R not clever enough to handle sum:
## diff by parms:
dpars <- c('k','s','a','b')
dD.fun <- function(theta){
  list2env(theta,envir = environment())
  -c(sum( (Pop*(exp(a+k+(s+b)*itz)-exp(k+ s*itz)))[itz>intT] ),
     sum( (Pop*itz*(exp(a+k+(s+b)*itz)-exp(k+s*itz)))[itz>intT] ),
     sum( (Pop*(exp(a+k+(s+b)*itz)))[itz>intT] ),
     sum( (Pop*itz*(exp(a+k+(s+b)*itz)))[itz>intT] ))
}

theta.woc <- list(
  intT = 8, #integer time before the post-intervention period (see graph above)
  itz=WT[!is.na(t),t],
  Pop=WT[!is.na(t),Pop],
  k=coef(mod.woc)['(Intercept)'],
  s=coef(mod.woc)['t'],
  a=coef(mod.woc)['Iafter'],
  b=coef(mod.woc)['t:Iafter']
)


## Check understanding of parameters:
## make data
tz <- seq(from=min(WT$year),to=max(WT$year),by=0.1) #time
tz <- seq(from=min(WT$t,na.rm=TRUE),to=max(WT$t,na.rm=TRUE),by=0.1) #time
It <- ifelse(tz>=1958,1,0)  #indicator
It <- ifelse(tz>=8,1,0)  #indicator
nap <- with(data=theta.woc,{exp(k + s*tz)}) #counterfactual
nap2 <- with(data=theta.woc,{exp(It*(a+b*tz))}) #ACF effect
napc <- (nap*nap2)[tz%%1==0]                    #with ACF at integer times
pop <- WT[!is.na(t),Pop]
TZ <- WT[,year]
TZ <- WT[,t]

## plot NOTE OK
png(file=here('exploratory/figures/delta.diagram.png'))
plot(WT[,.(t,tb)],ylim=c(0,3000))
lines(tz[tz%%1==0],napc*pop,col=2)
points(tz[tz%%1==0],nap[tz%%1==0]*pop,col=4)
points(WT[,t],
       exp(predict(mod.woc,newdata = WT)),
       col=2,pch='x')
dev.off()

## Checks: NOTE OK
blue <- nap[tz%%1==0]*pop
red <- napc*pop
sum(
  (blue-red)[TZ>8]
)
D.fun(theta.woc)              #our formula


## The names in the variance-covariance matrix are in the same order as the gradient terms, so:
Sig.woc <- vcov(mod.woc)
rownames(Sig.woc) <- colnames(Sig.woc) <- dpars

## Function to compute variance using delta approximation:
delta.se.woc <- function(theta){
  g <- dD.fun(theta)
  V <- t(g) %*% Sig.woc %*% (g)
  sqrt(V)
}
delta.se.woc(theta.woc) #test

## This would give the estimate as 95% CI
woc.result <- c(D.fun(theta.woc),
                D.fun(theta.woc) - 1.96*delta.se.woc(theta.woc),
                D.fun(theta.woc) + 1.96*delta.se.woc(theta.woc))
woc.result

cat(woc.result,file=here('exploratory/figures/notes.averted.afteronly.txt'))

## same as above but now including 'during'
prd <- exp(predict(mod.woc,newdata = WT))
WT
dfdur <- prd[8]-WT[t==8,tb]


## TODO include the up tick in notification and CI
