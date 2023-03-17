library(nimble)

species.names<-c("Americanwigeon","Canvasback","CinBlueTeal","Gadwall","Mallard","Northernpintail","Northernshoveler","Redhead","RuddyDuck")

for (sp in 1:length(species.names)){
 species.names<-c("Americanwigeon","Canvasback","CinBlueTeal","Gadwall","Mallard","Northernpintail","Northernshoveler","Redhead","RuddyDuck")
 load(paste0(species.names[sp],".RData"))

 source("ZINB.segment.nimble.R")

 duck<-nimbleModel(code = duckCode, name = "duck", constants = constants.nimble, data = data.nimble, inits = inits.nimble, calculate=FALSE)
 duckConf<-configureMCMC(duck, print = FALSE, thin=45)
 duckConf$addSampler(type = "RW_block",target = c("beta","mu.beta","gamma","mu.gamma"),control=list(adaptInterval = 50,adaptive=TRUE))
 duckConf$addSampler(type = "RW_block",target = c("r.us","theta"),control=list(adaptInterval = 50,adaptive=TRUE))
 duckConf$addMonitors(c('p','r.us','r','s2.p','N','theta','beta','mu.beta','gamma','mu.gamma','sum.lam','n','sumLogProb'))
 duckMCMC<- buildMCMC(duckConf)

 duckComp<-compileNimble(duck)
 duckComp.custom <- compileNimble(duckMCMC,project=duck,resetFunctions = TRUE)
 samplesList<-runMCMC(duckComp.custom, niter = 300000,  nburnin = 30000, nchains = 3, setSeed=c(10))

 save(constants.nimble,data.nimble,inits.nimble,samplesList,file=paste(species.names[sp],"early_LBC_Reproducing2.RData",sep="_"))
 keep(duck,duckComp,duckConf,duckMCMC,duckComp.custom,sp,sure=TRUE)
 gc()
}