library(ape)
library(magrittr)
library(lubridate)
library(skyspline)
library(treedater)
library(ggplot2)

## MCC tree from BEAST

mcc <- read.nexus("trees/Makona_1610_cds_ig.MCC.tree") 
mcc.sts <- mcc$tip.label %>% strsplit(.,"|",fixed=TRUE) %>% lapply(.,tail,1) %>% ymd %>% decimal_date
names(mcc.sts) <- mcc$tip.label
mcc.bdt <- DatedTree(mcc, mcc.sts, tol = 1/365)
mcc.t0 <- mcc.bdt$maxSampleTime - mcc.bdt$maxHeight
mcc.sky <- fit.skyspline.ml(mcc.bdt, death_rate_guess = 1/(15/365), t0=mcc.t0, y0_guess=1,  R0guess = 2, np_range = 2)

with(mcc.sky$demo.history, plot(times, pop.size, type='l'))
with(mcc.sky$demo.history, plot(times, reproduction.number, type='l'))

# cumulative infections 
GAMMA <- 1/(15/365); mcc.dx <- abs(diff( mcc.sky$demo.history$times)[1])
with(mcc.sky$demo.history, plot(rev(times), cumsum(mcc.dx*rev(pop.size)*rev(reproduction.number)*GAMMA), type='l'))
mcc.newinf <- with(mcc.sky$demo.history,mcc.dx*rev(pop.size)*rev(reproduction.number)*GAMMA)
mcc.maxinf <- with(mcc.sky$demo.history, date_decimal(times[mcc.newinf==max(mcc.newinf)]))
mcc.r <- with(mcc.sky$demo.history, max(reproduction.number))

## treedater

treefile <- "trees/Makona_1610_cds_ig_iqtree_allnni.treefile"
seqlen <- 18992
tre <- unroot(read.tree(treefile))
tipnames <- tre$tip.label
tipdates <- tipnames %>% strsplit(.,"_",fixed=TRUE) %>% lapply(.,tail,1) %>% ymd %>% decimal_date
names(tipdates) <- tipnames
sts <- tipdates[!is.na(tipdates)]
treed <- dater(tre, sts, s=seqlen, maxit = 100, minblen=1./365, abstol = .001, quiet = TRUE, temporalConstraints=TRUE, numStart=2, searchRoot=10, strictClock=FALSE, ncpu=10)

treed.bdt <- DatedTree(treed, sts, tol = 1.0/365)
treed.t0 <- treed.bdt$maxSampleTime - treed.bdt$maxHeight
treed.sky <- fit.skyspline.ml(treed.bdt, death_rate_guess = 1/(15/365), t0=treed.t0, y0_guess=1,  R0guess = 2, np_range = 2)
treed.dx <- abs(diff(treed.sky$demo.history$times)[1])
treed.newinf <- with(treed.sky$demo.history,treed.dx*rev(pop.size)*rev(reproduction.number)*GAMMA)
treed.maxinf <- with(treed.sky$demo.history, date_decimal(times[treed.newinf==max(treed.newinf)]))
treed.r <- with(treed.sky$demo.history, max(reproduction.number))

with(treed.sky$demo.history, plot(times, pop.size, type='l'))
with(treed.sky$demo.history, plot(times, reproduction.number, type='l'))


### Now add mutations to all edges

tre2 <- tre
set.seed(1234)
tre2$edge.length[1:length(tre2$edge.length)] <- tre2$edge.length[1:length(tre2$edge.length)]+runif(length(tre2$edge.length),min=0,max=1./seqlen)
treed.2 <- dater(tre2, sts, s=seqlen, maxit = 100, minblen=1./365, abstol = .001, quiet = TRUE, temporalConstraints=TRUE, numStart=2, searchRoot=10, strictClock=FALSE, ncpu=10)

treed.2.bdt <- DatedTree(treed.2, sts, tol = 1.0/365)
treed.2.t0 <- treed.2.bdt$maxSampleTime - treed.2.bdt$maxHeight
treed.2.sky <- fit.skyspline.ml(treed.2.bdt, death_rate_guess = 1/(15/365), t0=treed.2.t0, y0_guess=1,  R0guess = 2, np_range = 2)
treed.2.dx <- abs(diff( treed.2.sky$demo.history$times)[1])
treed.2.newinf <- with(treed.2.sky$demo.history, treed.2.dx*rev(pop.size)*rev(reproduction.number)*GAMMA)
treed.2.maxinf <- with(treed.2.sky$demo.history, date_decimal(times[treed.2.newinf==max(treed.2.newinf)]))
treed.2.r <- with(treed.2.sky$demo.history, max(reproduction.number))

### Mutations to tips only

tre3 <- tre
set.seed(1234)
tre3$edge.length[1:length(tre3$tip.label)] <- tre3$edge.length[1:length(tre3$tip.label)]+runif(length(tre3$tip.label),min=0,max=1./seqlen)
treed.3 <- dater(tre3, sts, s=seqlen, maxit = 100, minblen=1./365, abstol = .001, quiet = TRUE, temporalConstraints=TRUE, numStart=2, searchRoot=10, strictClock=FALSE, ncpu=10)
treed.3.bdt <- DatedTree(treed.3, sts, tol = 1.0/365)
treed.3.t0 <- treed.3.bdt$maxSampleTime - treed.3.bdt$maxHeight
treed.3.sky <- fit.skyspline.ml(treed.3.bdt, death_rate_guess = 1/(15/365), t0=treed.3.t0, y0_guess=1,  R0guess = 2, np_range = 2)
treed.3.dx <- abs(diff( treed.3.sky$demo.history$times)[1])
treed.3.newinf <- with(treed.3.sky$demo.history, treed.3.dx*rev(pop.size)*rev(reproduction.number)*GAMMA)
treed.3.maxinf <- with(treed.3.sky$demo.history, date_decimal(times[treed.3.newinf==max(treed.3.newinf)]))
treed.3.r <- with(treed.3.sky$demo.history, max(reproduction.number))

## node.dating

tre.rtt <- rtt(tre,tipdates,objective="rsquared")
td <- tre.rtt$tip.label %>% strsplit(.,"_",fixed=TRUE) %>% lapply(.,tail,1) %>% ymd %>% decimal_date
tre.rtt2 <- tre.rtt
#set.seed(1235)
#tre.rtt2$edge.length[1:length(tre.rtt2$edge.length)] <- tre.rtt2$edge.length[1:length(tre.rtt2$edge.length)]+runif(length(tre.rtt2$edge.length),min=0,max=1./seqlen)
mu <- estimate.mu(tre.rtt2, td)
tre.rtt.dates <- estimate.dates(tre.rtt2, td, mu, nsteps = 100)
nd <- tre.rtt2
nd$edge.length <- tre.rtt.dates[tre.rtt2$edge[, 2]] - tre.rtt.dates[tre.rtt2$edge[, 1]]

nd.sts <- nd$tip.label %>% strsplit(.,"_",fixed=TRUE) %>% lapply(.,tail,1) %>% ymd %>% decimal_date
names(nd.sts) <- nd$tip.label
nd.bdt <- DatedTree(nd, nd.sts, tol = 1/365)
nd.t0 <- nd.bdt$maxSampleTime - nd.bdt$maxHeight
nd.sky <- fit.skyspline.ml(nd.bdt, death_rate_guess = 1/(15/365), t0=nd.t0, y0_guess=1,  R0guess = 2, np_range = 2)
nd.r <- with(nd.sky$demo.history, max(reproduction.number))

with(nd.sky$demo.history, plot(times, pop.size, type='l'))
with(nd.sky$demo.history, plot(times, reproduction.number, type='l'))

# cumulative infections 
GAMMA <- 1/(15/365); nd.dx <- abs(diff( nd.sky$demo.history$times)[1])
with(nd.sky$demo.history, plot(rev(times), cumsum(nd.dx*rev(pop.size)*rev(reproduction.number)*GAMMA), type='l'))
tail(with(nd.sky$demo.history, cumsum(nd.dx*rev(pop.size)*rev(reproduction.number)*GAMMA)),1)
nd.dx <- abs(diff(nd.sky$demo.history$times)[1])
nd.newinf <- with(nd.sky$demo.history,nd.dx*rev(pop.size)*rev(reproduction.number)*GAMMA)
nd.maxinf <- with(nd.sky$demo.history, date_decimal(times[nd.newinf==max(nd.newinf)]))


## lsd

lsd <- read.tree("results/Makona_1610_cds_ig.lsd.date.newick") 
lsd.sts <- lsd$tip.label %>% strsplit(.,"_",fixed=TRUE) %>% lapply(.,tail,1) %>% ymd %>% decimal_date
names(lsd.sts) <- lsd$tip.label
lsd.bdt <- DatedTree(lsd, lsd.sts, tol = 1/365)
lsd.t0 <- lsd.bdt$maxSampleTime - lsd.bdt$maxHeight
lsd.sky <- fit.skyspline.ml(lsd.bdt, death_rate_guess = 1/(15/365), t0=lsd.t0, y0_guess=1,  R0guess = 2, np_range = 2)
lsd.r <- with(lsd.sky$demo.history, max(reproduction.number))

with(lsd.sky$demo.history, plot(times, pop.size, type='l'))
with(lsd.sky$demo.history, plot(times, reproduction.number, type='l'))

# you can also get cumulative infections using this weird incantation: 
GAMMA <- 1/(15/365); lsd.dx <- abs(diff( lsd.sky$demo.history$times)[1])
with(lsd.sky$demo.history, plot(rev(times), cumsum(lsd.dx*rev(pop.size)*rev(reproduction.number)*GAMMA), type='l'))
tail(with(lsd.sky$demo.history, cumsum(lsd.dx*rev(pop.size)*rev(reproduction.number)*GAMMA)),1)

lsd.newinf <- with(lsd.sky$demo.history, lsd.dx*rev(pop.size)*rev(reproduction.number)*GAMMA)
lsd.maxinf <- with(lsd.sky$demo.history, date_decimal(times[lsd.newinf==max(lsd.newinf)]))

## Cumulative infections

lsd.cs <- with(lsd.sky$demo.history,
     data.frame(Time=rev(times),Cumulative=cumsum(lsd.dx*rev(pop.size)*rev(reproduction.number)*GAMMA)))
mcc.cs <- with(mcc.sky$demo.history,
               data.frame(Time=rev(times),Cumulative=cumsum(mcc.dx*rev(pop.size)*rev(reproduction.number)*GAMMA)))
treed.cs <- with(treed.sky$demo.history,
                    data.frame(Time=rev(times),Cumulative=cumsum(treed.dx*rev(pop.size)*rev(reproduction.number)*GAMMA)))
treed.2.cs <- with(treed.2.sky$demo.history,
                   data.frame(Time=rev(times),Cumulative=cumsum(treed.2.dx*rev(pop.size)*rev(reproduction.number)*GAMMA)))
treed.3.cs <- with(treed.3.sky$demo.history,
                    data.frame(Time=rev(times),Cumulative=cumsum(treed.3.dx*rev(pop.size)*rev(reproduction.number)*GAMMA)))
nd.cs <- with(nd.sky$demo.history,
                    data.frame(Time=rev(times),Cumulative=cumsum(nd.dx*rev(pop.size)*rev(reproduction.number)*GAMMA)))

obs.cs <- read.table("ebov_cs.txt",header=T,row.names=NULL)
obs.cs$Time <- decimal_date(dmy(obs.cs$Date))
obs.cs$Cumulative <- obs.cs$Cases
obs.cs <- obs.cs[obs.cs$Time<max(lsd.cs$Time),]

p <- ggplot()+geom_line(data=treed.cs,aes(x=Time,y=Cumulative,colour="treedater+ML tree"))+
  geom_line(data=treed.2.cs,aes(x=Time,y=Cumulative,colour="treedater+edge mutations"))+
  geom_line(data=treed.3.cs,aes(x=Time,y=Cumulative,colour="treedater+tip mutations"))+
  geom_line(data=mcc.cs,aes(x=Time,y=Cumulative,colour="BRMC"))+
  geom_line(data=lsd.cs,aes(x=Time,y=Cumulative,colour="QPD"))+
  geom_line(data=nd.cs,aes(x=Time,y=Cumulative,colour="node dating"))+
  geom_line(data=obs.cs,aes(x=Time,y=Cumulative,colour="Observed"))+
  xlab("Time")+ylab("Cumulative number of infections")+
  scale_colour_manual(name="Model",
                      values=c("treedater+ML tree"="red", "treedater+edge mutations"="orange","treedater+tip mutations"="blue","BRMC"="purple", "QPD"="green", "node dating"="cyan","Observed"="black"))
ggsave("cumulative.pdf",p,width=8,height=5,units="in")

## Supplemental figure of edge lengths

el <- data.frame(EdgeLength=c(mcc$edge.length,treed$edge.length),Tree=c(rep("BRMC MCC",length(mcc$edge.length)),rep("treedater+ML tree",length(treed$edge.length))))
p2 <- ggplot(el,aes(x=EdgeLength))+geom_histogram(aes(y=..density..),binwidth=0.01)+facet_wrap(~Tree)+xlab("Edge length (years)")+ylab("Density")
ggsave("edgelengths.pdf",p2,width=8,height=5,units="in")

## LTT

ltt <- function(tre,mrsd){
  lineages <- ltt.plot.coords(tre,backward=TRUE)
  lineages <- as.data.frame(lineages)
  names(lineages) <- c("Time","Lineages")
  lineages$Time <- mrsd+lineages$Time
  lineages
}

treed.ltt <- ltt(treed.bdt,decimal_date(ymd("2015-10-24")))
mcc.ltt <- ltt(mcc.bdt,decimal_date(ymd("2015-10-24")))
p3 <- ggplot()+geom_line(data=treed.ltt,aes(x=Time,y=Lineages,color="treedater"))+geom_line(data=mcc.ltt,aes(x=Time,y=Lineages,color="BRMC"))+scale_colour_manual(name="Model",
                                                                                                                                                                       values=c("treedater"="red","BRMC"="purple"))
ggsave("ltt.pdf",p3,width=8,height=5,units="in")
