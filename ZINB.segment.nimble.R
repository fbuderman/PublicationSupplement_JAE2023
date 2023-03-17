duckCode <- nimbleCode(
{
 sumLogProb ~ dnorm(0,1) 

 r.us~dnorm(r.prior,1000)

 N~dunif(0,20)

 p~dbeta(1,1)

 theta~dunif(0,1)

 for (l in 1:X.cov.num){
  tau.beta[l]~dgamma(0.0001,0.0001)
  mu.beta[l]~dnorm(0,0.001)
 }

 for (k in 1:Z.cov.num){
  tau.gamma[k]~dgamma(0.0001,0.0001)
  mu.gamma[k]~dnorm(0,0.001)
 }

 for (s in 1:S){
  tau.p[s]~dgamma(0.001,0.001)
  s2.p[s]<-1/tau.p[s]

  log.n.0[s]~dnorm(0,0.1)
  n.0[s]<-exp(log.n.0[s])

  for (l in 1:X.cov.num){
   beta[s,l]~dnorm(mu.beta[l],tau.beta[l])
  } 

  for (k in 1:Z.cov.num){
   gamma[s,k]~dnorm(mu.gamma[k],tau.gamma[k])
  } 
 }

 for (i in 1:obs){
  y[i] ~ dnegbin(y.mean[i],N)
  y.mean[i]<-N/(N+zero[i]*lam.i[i])-1e-10*(1-zero[i])
  zero[i] ~ dbern(p)
  log.lam.i[i]<-log(n[seg.year[i],strata.seg[i]]/I.in.S[i])+inprod(X.segment[i,1:X.cov.num],beta[strata.seg[i],1:X.cov.num])
  lam.i[i]<-exp(log.lam.i[i])
 }

 for (m in 1:SY){
  sum.lam[m]<-sum(lam.i[start.end.strata.segment[1,m]:start.end.strata.segment[2,m]])
  log.sum.lam[m]<-log(sum(lam.i[start.end.strata.segment[1,m]:start.end.strata.segment[2,m]]))
 }

 for (s in 1:S){
  r[s]<-r.us+sum(Z.attributes[1:Z.cov.num]*gamma[s,1:Z.cov.num])
  log.n[1,s]~dnorm(log.n.mean[1,s],tau.p[s])
  log.n.mean[1,s]<-r[s]/(1-theta)+(theta*log.n.0[s])+inprod(Z[1,1:(Z.cov.num-1),s],gamma[s,1:(Z.cov.num-1)])
  n[1,s]<-exp(log.n[1,s])
  for (t in 2:T){
   log.n[t,s]~dnorm(log.n.mean[t,s],tau.p[s])
   log.n.mean[t,s]<-r[s]/(1-theta)+(theta*log.sum.lam[id[t-1,s]])+inprod(Z[t,1:Z.cov.num,s],gamma[s,1:Z.cov.num])
   n[t,s]<-exp(log.n[t,s])
  }
 }
}
)