#prepare environment
library(jagsUI)
library(MASS)
setwd("~/Nest Survival Simulation/simulations")

#DEFINE SIMULATION FUNCTION####
nest.fate.fn.fixedeffs<-function(n.marked = 200, beta.a.ex=1.8, alpha.a=-8.3, n.ex=0.5, JAGSmodel="ToolModelWidePriors.txt"){

#set nest initiation dates
init<-c(1:n.marked)
for(i in 1:length(init)){
	init[i]<-rnbinom(1,84,mu=32)  #from nest initiation data, standardized to first nest of all data = day 1
	}

n.occasions<-120  #length of breeding season

#intercept and exclosure effects on predation (p) and abandonment (a) and flooding (f)  #should I also add two covariates for flooding? test = 2 covariates per fate
alpha.p<- -4.2
alpha.f<- -6.4
beta.p.ex<- -2.1

#sample nests to be exclosed
exclose.nest<-rep(rbinom(n.marked,1,n.ex)) 
ex<-matrix(0,ncol=n.occasions, nrow=n.marked) #initialize empty exclosure matrix

#fill in exclosure matrix - nests must be at least 4 days old, plus additional variation in timing of placement (1-4 extra days)
for(i in 1:n.marked){
	ex[i,(init[i]+4+sample(1:4,1)):n.occasions]<-exclose.nest[i]
}

#prepare fate matrices and linear predictors
pred<-matrix(0,ncol=n.occasions, nrow=n.marked)
aban<-matrix(0,ncol=n.occasions, nrow=n.marked)
flood<-matrix(0,ncol=n.occasions, nrow=n.marked)
surv<-matrix(0,ncol=n.occasions, nrow=n.marked)

ctp<-exp(alpha.p+beta.p.ex*ex)
cta<-exp(alpha.a+beta.a.ex*ex)
ctf<-exp(alpha.f)
cts<-1

den<-ctp+cta+cts+ctf
survp<-cts/den

pred<-((ctp/den)/(1-survp)*(1-survp))  
aban<-((cta/den)/(1-survp)*(1-survp))
flood<-((ctf/den)/(1-survp)*(1-survp))
surv<-survp  #daily survival probability

PHI<-array(NA, dim=c(n.marked,n.occasions,4))
PHI[,,1]<-surv
PHI[,,2]<-pred
PHI[,,3]<-flood
PHI[,,4]<-aban
hatch<-34

Fate.mat<-matrix(0, ncol=n.occasions, nrow=n.marked)   #empty matrix of nest fates
final.fate<-rep(0,n.marked)

#Fill in fate matrix
for(i in 1:n.marked){
  Fate.mat[i,init[i]]<-1  #put a 1 at release occasion (nest initiation date)
  if (init[i]==n.occasions) next
  for(t in (init[i]+1):n.occasions){
	status<-which(rmultinom(1,1,PHI[i,t-1,])==1) 
	Fate.mat[i,t]<-status 
	if (status==2) break #move on if predation
	if (status==3) break  #flood
	if (status==4) break #abandonment
	if (t-(init[i]-1)==hatch) break  #allow nests to hatch
	}#t
	final.fate[i]<-max(Fate.mat[i,])
 }#i

###observation process: probability of discovering a nest, and checking less frequently than every day###
disc=0.9  #probability of discovering an active nest
freq=1  #monitor every day for simplicity 
p=disc*freq  #average daily probability of discovering a nest
encounter.history<-matrix(0,ncol=n.occasions, nrow=n.marked)
discovery.date<-rep(NA, n.marked)  
last.check<-rep(NA, n.marked)

#discover nests
for(i in 1:n.marked){
	for(t in (init[i]):n.occasions){
		if (Fate.mat[i,t]==1) encounter.history[i,t]=rbinom(1,1,p)
		if (encounter.history[i,t]==1) break
		}#t 
		if (sum(encounter.history[i,])>0) discovery.date[i]<-which(encounter.history[i,]==1)
		else discovery.date[i]<-n.occasions #problems created if assigned NA
}#i

#fill in remaining encounter history
for(i in 1:n.marked){
	for(t in (discovery.date[i]):n.occasions){
		if (discovery.date[i]==n.occasions) break
		if (Fate.mat[i,t]==1) encounter.history[i,t]=rbinom(1,1,freq)
	}#t
        if (sum(encounter.history[i,])==0) {final.fate[i]<-0
	    last.check[i] <- 1} #deal with nests never encountered to avoid warnings
	else 
	{last.check[i]<-(max(which(encounter.history[i,]>0))+sample(0:2,1)) }  #record final nest fate 0-2 days after it happens	
	if (last.check[i] > n.occasions) last.check[i]<-n.occasions
	encounter.history[i,last.check[i]]<-final.fate[i]
      
}#i

#rearrange data into vectorized format with covariates
fate.vec<-c(t(encounter.history))
exclose.vec<-c(t(ex))
days<-rep(1:n.occasions,n.marked)
nestID<-rep(c(1:n.marked), each=n.occasions)
nestdata<-as.data.frame(cbind(fate.vec,days,nestID,exclose.vec))

nestdata=nestdata[nestdata$fate.vec!=0,]  #select only dates with nest checks

#calculate exposure days
expose<-rep(NA, length(nestdata$nestID))
for(i in 2:length(nestdata$nestID)){
  if (nestdata$nestID[i]==nestdata$nestID[i-1]) {expose[i]=(nestdata$days[i]-nestdata$days[i-1])}
  else {expose[i]=NA}
}

#bind to data and remove first encounters 
nestdata$expose<-expose
nestdata<-nestdata[!(is.na(nestdata$expose)),]

#prepare fate matrix for JAGS
n<-length(nestdata$expose)
Surv=rep(0,n)
Aban=rep(0,n)
Pred=rep(0,n)
Flood=rep(0,n)

for (i in 1:n){
Surv[i][nestdata$fate.vec[i]==1]=1
Aban[i][nestdata$fate.vec[i]==4]=1
Pred[i][nestdata$fate.vec[i]==2]=1
Flood[i][nestdata$fate.vec[i]==3]=1
}

Fate<-cbind(Surv, Aban, Pred, Flood)

jags.data<-list(ex=nestdata$exclose.vec,n=n, interval=nestdata$expose, Fate=Fate)

#set initial values
inits <- function() {list(
		    alpha.p = rnorm(1, 0, 1), alpha.a = rnorm(1, 0, 1), alpha.f=rnorm(1,0,1),
			  beta.a.ex = rnorm(1, 0, 1), beta.p.ex = rnorm(1,0,1))}

parameters <- c("alpha.p","alpha.a","alpha.f","beta.a.ex","beta.p.ex")

out <- autojags(jags.data, inits, parameters, model.file=JAGSmodel, n.chains=3, parallel=TRUE)

vcovAban <- cov(cbind(out$sims.list$alpha.a, out$sims.list$beta.a.ex))
vcovPred <- cov(cbind(out$sims.list$alpha.p, out$sims.list$beta.p.ex))

max.Rhat<-max(as.numeric(lapply(out$Rhat, function(x)x[which.max(x)])))  #check for successful convergence of all parameters
converge.success <- T
if (max.Rhat>1.1) {converge.success<-F} 

output<-list(out=out, aban=sum(Aban), pred=sum(Pred), flood=sum(Flood), vcovAban=vcovAban, vcovPred=vcovPred,
             converge.success=converge.success)
return(output)
} #nest fate simulation function

#HATCH CALC FXNS####
h.ex.calc<-function(parms){
  
  lin.f.ex<-exp(parms$alpha.f)
  lin.a.ex<-exp(parms$alpha.a+ parms$beta.a.ex)
  lin.p.ex<-exp(parms$alpha.p  + parms$beta.p.ex)
  survp.ex<-1/(lin.f.ex+lin.a.ex+lin.p.ex+1)
  pred.p.ex<-((lin.p.ex/(lin.f.ex+lin.a.ex+lin.p.ex+1))/(1-survp.ex)*(1-survp.ex^34))
  a<-((lin.a.ex/(lin.f.ex+lin.a.ex+lin.p.ex+1))/(1-survp.ex)*(1-survp.ex^34))
  flood.p.ex<-((lin.f.ex/(lin.f.ex+lin.a.ex+lin.p.ex+1))/(1-survp.ex)*(1-survp.ex^34))
  
  h<-survp.ex^34   #probability of period survival exclosed nest
  o<-flood.p.ex+pred.p.ex 
  m=parms$m
  r2=parms$r2
  r3=parms$r3
  h.ex<-h + (o*r2*h + a*(1-m)*r2*h) + o*r2*(o*r3*h + a*(1-m)*r3*h) + a*(1-m)*r2*(o*r3*h+a*(1-m)*r3*h)
  output<-list(h.ex=h.ex, a=a)
  return(output)
}

h.un.calc<-function(parms){
  lin.f<-exp(parms$alpha.f)
  lin.a<-exp(parms$alpha.a)
  lin.p<-exp(parms$alpha.p)
  survp<-1/(lin.f+lin.a+lin.p+1)   #daily survival probability
  
  #period fate probabilities for each nest
  pred.p<-((lin.p/(lin.f+lin.a+lin.p+1))/(1-survp)*(1-survp^34))  #predation unexclosed
  a<-((lin.a/(lin.f+lin.a+lin.p+1))/(1-survp)*(1-survp^34))  #abandonment unexclosed
  flood.p<-((lin.f/(lin.f+lin.a+lin.p+1))/(1-survp)*(1-survp^34))  #flooding probability
  o<-flood.p+pred.p
  h<-survp^34  #probability of period survival unexclosed nest
  m=parms$m
  r2=parms$r2
  r3=parms$r3
  h.un <- h + (o*r2*h + a*(1-m)*r2*h) + o*r2*(o*r3*h + a*(1-m)*r3*h) + a*(1-m)*r2*(o*r3*h+a*(1-m)*r3*h) #revised; mistake somewhere in original
  output<-list(h.un=h.un, a=a)
  return(output)
}
#logodds####
logodds=function(x){log(x/(1-x))}

#LAMBDA CALC FXN####
lambda.calc <- function(parms, sd.parms, eta.a=0, eta.p=0, n.ex=0, n.iter=1000){
  lambda <- rep(NA, n.iter)
  #run iterations
  for (i in 1:n.iter){
    #draw stochastic values
    parms.iter <- list(
      survmean=plogis(mvrnorm(1,c(qlogis(parms$Phij),qlogis(parms$Phia)),sd.parms$vcovSurv)),
      f = plogis(rnorm(1,qlogis(parms$f),sd.parms$sd.f)),
      ys = rnorm(1,logodds(parms$ys),sd.parms$ys), #probability of second-year bird nesting
      yt = parms$yt, #probability of third-year bird nesting fixed at 0.99
      r2 = rnorm(1, logodds(parms$r2), sd.parms$r2),
      r3 = rnorm(1, logodds(parms$r3), sd.parms$r3),
      m = plogis(rnorm(1,qlogis(parms$m), sd.parms$m)),
      alpha.f = rnorm(1, parms$alpha.f, sd.parms$alpha.f),
      aban.coeff = mvrnorm(1, c(parms$alpha.a, parms$beta.a.ex), sd.parms$vcovAban),
      pred.coeff = mvrnorm(1, c(parms$alpha.p, parms$beta.p.ex), sd.parms$vcovPred)
    )
    parms.iter$alpha.a <- parms.iter$aban.coeff[1] + eta.a
    parms.iter$beta.a.ex <- parms.iter$aban.coeff[2]
    parms.iter$alpha.p <- parms.iter$pred.coeff[1] + eta.p
    parms.iter$beta.p.ex <- parms.iter$pred.coeff[2]
    
    #baseline abandonment rate (unexclosed)
    a.p.un<-h.un.calc(parms=parms.iter)$a
    
    #hatch and exclosure-related abandonment probabilities
    h.un<-h.un.calc(parms=parms.iter)$h.un
    h.ex<-h.ex.calc(parms=parms.iter)$h.ex
    a.p.ex<-h.ex.calc(parms=parms.iter)$a
    
    #fecundity
    Fa<-parms.iter$yt*(2*parms$E*(n.ex*h.ex + (1-n.ex)*h.un))
    Fs<-parms.iter$ys*(2*parms$E*(n.ex*h.ex + (1-n.ex)*h.un))
    
    #breeding-season mortality with abandonment
    m<-parms.iter$m; r2<-parms.iter$r2; r3<-parms.iter$r3
    m.a.un=parms.iter$yt*(a.p.un*m+a.p.un*(1-m)*r2*(a.p.un*m+a.p.un*(1-m)*r3*a.p.un*m)) 
    m.s.un=parms.iter$ys*(a.p.un*m+a.p.un*(1-m)*r2*(a.p.un*m+a.p.un*(1-m)*r3*a.p.un*m))
    m.a.ex=parms.iter$yt*(a.p.ex*m+a.p.ex*(1-m)*r2*(a.p.ex*m+a.p.ex*(1-m)*r3*a.p.ex*m))  
    m.s.ex=parms.iter$ys*(a.p.ex*m+a.p.ex*(1-m)*r2*(a.p.ex*m+a.p.ex*(1-m)*r3*a.p.ex*m))

    #annual survival rates	
    Phij.w <- parms.iter$survmean[1]^(10/12) #winter post-fledging juvenile survival
    Phia.w <- parms.iter$survmean[2]^(10/12) #winter adult survival
    Phij.b <- parms.iter$survmean[1]^(2/12)  #second-year breeding season survival
    Phia.b <- parms.iter$survmean[2]^(2/12)  #ASY breeding survival
    Phij.b <- Phij.b/(1-m.s.un)*(1-m.s.ex*n.ex)  #add in probability of surviving exclosure-related abandonment
    Phia.b <- Phia.b/(1-m.a.un)*(1-m.a.ex*n.ex)
    Phij.ann <- Phij.b*Phij.w 
    Phia.ann <- Phia.b*Phia.w
    
    #matrix calculations
    f<-parms.iter$f
    Lesmat<-matrix(c(Phij.ann*Fs*f,Fa*Phia.ann,Phij.ann*f,Phia.ann),2,2,byrow=TRUE)  #found mistake in original; entry [2,1] had adult survival, not juv
    lambda[i]<-eigen(Lesmat)$values[1]
  } #n.iter
  output <- list(lambda=lambda)
  return(output)
}

#MEAN PARAMETERS####

mean.parms <- list(
  yt=0.99,  #probability of third-year bird nesting
  ys=0.68,  #probability of second-year bird nesting
  Phij=0.52,  
  Phia=0.74,
  f=0.4,
  m=0.34,
  E = 0.94,
  r2 = 0.7,
  r3 = 0.7,
  beta.a.ex = 1.284,
  beta.p.ex = -2.532,
  alpha.a = -7.447,
  alpha.p = -4.257,
  alpha.f = -6.040
)	

mean.alpha.a = -7.447
mean.alpha.p = -4.257
CV=0.1 #coefficient of variation for most stochastic parameters

#PARAMETER VARIANCES####	
sd.parms.init <- list(
  vcovSurv=matrix(c((logodds(mean.parms$Phij)*CV)^2,logodds(mean.parms$Phij)*CV*logodds(mean.parms$Phia)*CV, 
                    logodds(mean.parms$Phij)*CV*logodds(mean.parms$Phia)*CV, (logodds(mean.parms$Phia)*CV))^2,2,2),  #had sd in here earlier, not variance
  sd.f = abs(logodds(mean.parms$f))*0.06,
  yt = logodds(mean.parms$yt)*CV,
  ys = logodds(mean.parms$ys)*CV,
  r2 = abs(logodds(mean.parms$r2))*CV,
  r3 = abs(logodds(mean.parms$r3))*CV,
  m = abs(logodds(mean.parms$m))*1.04,  
  beta.a.ex = 0.348,
  beta.p.ex = 0.272,
  alpha.a = 0.386,
  alpha.p = 0.116,
  alpha.f = 0.222,
  vcovAban = matrix(c(0.386^2, -0.092, -0.092, 0.348^2),2,2),
  vcovPred = matrix(c(0.116^2, -0.007, -0.007, 0.272^2),2,2)
)

#SIMULATIONS####
#initiate vectors to store output results
n=100
sample.sizes <- c(5,10,20,50,100)
n.samples <- length(sample.sizes)
alpha.a.vec <- matrix(NA, nrow=n, ncol=n.samples)
alpha.p.vec <- matrix(NA, nrow=n, ncol=n.samples)
alpha.f.vec <- matrix(NA, nrow=n, ncol=n.samples)
beta.a.ex.vec <- matrix(NA, nrow=n, ncol=n.samples)
beta.p.ex.vec <- matrix(NA, nrow=n, ncol=n.samples)
aban.count <- matrix(NA, nrow=n, ncol=n.samples)
pred.count <- matrix(NA, nrow=n, ncol=n.samples)
flood.count <- matrix(NA, nrow=n, ncol=n.samples)
converge <- matrix(NA, nrow=n, ncol=n.samples)
lambda.ex<- array(NA, dim=c(n,1000,n.samples))
lambda.un<-array(NA, dim=c(n,1000,n.samples))

#scenario 1: abandonment greater than average####
alpha.a = mean.parms$alpha.a + 2
beta.a.ex = mean.parms$beta.a.ex + 0.2*mean.parms$beta.a.ex

for (i in 1:n){ 
  for(s in 1:n.samples){
	sim <- nest.fate.fn.fixedeffs(alpha.a=alpha.a, beta.a.ex=beta.a.ex, n.marked=sample.sizes[s], JAGSmodel = "ToolModelInformPriors.txt")  
	alpha.a.vec[i,s] <- sim$out$mean$alpha.a
	alpha.p.vec[i,s] <- sim$out$mean$alpha.p
	alpha.f.vec[i,s] <- sim$out$mean$alpha.f
	beta.p.ex.vec[i,s] <-sim$out$mean$beta.p.ex
	beta.a.ex.vec[i,s] <- sim$out$mean$beta.a.ex
	aban.count[i,s] <- sim$aban
	pred.count[i,s] <- sim$pred
	flood.count[i,s] <- sim$flood
	converge[i,s]<- sim$converge.success
	mean.parms$beta.a.ex=sim$out$mean$beta.a.ex
	mean.parms$beta.p.ex=sim$out$mean$beta.p.ex
	mean.parms$alpha.a=sim$out$mean$alpha.a
	mean.parms$alpha.p=sim$out$mean$alpha.p
	mean.parms$alpha.f=sim$out$mean$alpha.f
	sd.parms.init$vcovAban=sim$vcovAban
	sd.parms.init$vcovPred=sim$vcovPred
	sd.parms.init$alpha.f=sim$out$sd$alpha.f
	sd.parms.init$alpha.a=sim$out$sd$alpha.a
	sd.parms.init$alpha.p=sim$out$sd$alpha.p
  lambda.ex[i,,s] <- lambda.calc(parms=mean.parms, sd.parms=sd.parms.init, n.ex=1)$lambda
  lambda.un[i,,s] <- lambda.calc(parms=mean.parms, sd.parms=sd.parms.init, n.ex=0)$lambda
	cat("Finished",i,"of",n,"runs in ", s, "of", n.samples, "simulations") 
  }
}

#decision criteria
rapid.growth <- array(NA, dim=c(n,2,n.samples))
growth <- array(NA, dim=c(n,2,n.samples))
decline <- array(NA, dim=c(n,2,n.samples))
rapid.decline <- array(NA, dim=c(n,2,n.samples))

for(i in 1:n){
  for(s in 1:n.samples){
    rapid.growth[i,1,s] <- mean(lambda.ex[i,,s]>1.05)
    rapid.growth[i,2,s] <- mean(lambda.un[i,,s]>1.05)
    growth[i,1,s] <- mean(lambda.ex[i,,s]>1)
    growth[i,2,s] <- mean(lambda.un[i,,s]>1)
    decline[i,1,s] <- mean(lambda.ex[i,,s]<1)
    decline[i,2,s] <- mean(lambda.un[i,,s]<1)
    rapid.decline[i,1,s] <- mean(lambda.ex[i,,s]<1.05) #column 1 = exclosed
    rapid.decline[i,2,s] <- mean(lambda.un[i,,s]<1.05) #column 2 = not exclosed
  }
}

#"correct" decision based on 100 nests (re-run with 1000?) is not exclosed; 0.99 probability of rapid decline with ex
mean(rapid.decline[,1,5])
#[1] 0.99003
#vs. 0.51 probability of rapid decline without exclosures:
mean(rapid.decline[,2,5])
#[1] 0.51671

#create decision matrices; compare probabilities of growth and decline with and w/out exclosures
#decision is based on greater growth probability and smaller decline probability
decide.growth <- decide.rapid.decline <- decide.decline <- decide.rapid.growth <- matrix(NA, nrow=n, ncol=n.samples)
decide.rapid.decline[,]<- rapid.decline[,1,]>rapid.decline[,2,]
decide.decline[,]<-decline[,1,]>decline[,2,]
decide.growth[,]<-growth[,1,]<growth[,2,]
decide.rapid.growth[,]<-rapid.growth[,1,]<rapid.growth[,2,]

#probability of making correct decision with n nests based on probability of rapid decline
mean(decide.rapid.decline[,1]) #5 nests
#[1] 0.55

#compare to a single reference run with 1000 nests - "true decision" (ran but did not preserve all code)
n=1 
ref.alpha.a.vec <- rep(NA, n)
ref.alpha.p.vec <- rep(NA, n)
ref.alpha.f.vec <- rep(NA, n)
ref.beta.a.ex.vec <- rep(NA, n)
ref.beta.p.ex.vec <- rep(NA, n)
ref.aban.count <- rep(NA, n)
ref.pred.count <- rep(NA, n)
ref.flood.count <- rep(NA, n)
ref.converge <- rep(NA, n)
ref.lambda.ex<- matrix(NA, nrow=n, ncol=1000)
ref.lambda.un<-matrix(NA, nrow=n, ncol=1000)



#SIMULATIONS WITH WIDE PRIORS####
ref.wide.lambda.ex <-matrix(NA, nrow=n, ncol=1000)
ref.wide.lambda.un <-matrix(NA, nrow=n, ncol=1000)

for (i in 1:n){ 
  sim.wide <- nest.fate.fn.fixedeffs(n.marked=1000, n.ex=0.7,alpha.a=alpha.a, beta.a.ex=beta.a.ex, JAGSmodel = "ToolModelWidePriors.txt")  
  mean.parms$beta.a.ex=sim.wide$out$mean$beta.a.ex
  mean.parms$beta.p.ex=sim.wide$out$mean$beta.p.ex
  mean.parms$alpha.a=sim.wide$out$mean$alpha.a
  mean.parms$alpha.p=sim.wide$out$mean$alpha.p
  mean.parms$alpha.f=sim.wide$out$mean$alpha.f
  sd.parms.init$vcovAban=sim.wide$vcovAban
  sd.parms.init$vcovPred=sim.wide$vcovPred
  sd.parms.init$alpha.f=sim.wide$out$sd$alpha.f
  sd.parms.init$alpha.a=sim.wide$out$sd$alpha.a
  sd.parms.init$alpha.p=sim.wide$out$sd$alpha.p
  ref.wide.lambda.ex[i,] <- lambda.calc(parms=mean.parms, sd.parms=sd.parms.init, n.ex=1)$lambda
  ref.wide.lambda.un[i,] <- lambda.calc(parms=mean.parms, sd.parms=sd.parms.init, n.ex=0)$lambda
  cat("Finished",i,"of",n,"runs")  
}


#did not finish running - extreme estimates causing failure in eigen(lambda)
wide.converge <- matrix(NA, nrow=n, ncol=n.samples)
wide.lambda.ex<- array(NA, dim=c(n,1000,n.samples))
wide.lambda.un<-array(NA, dim=c(n,1000,n.samples))

for (i in 1:n){ 
  for(s in 1:n.samples){
    sim <- nest.fate.fn.fixedeffs(alpha.a=alpha.a, beta.a.ex=beta.a.ex, n.marked=sample.sizes[s], 
                                  JAGSmodel = "ToolModelWidePriors.txt")  
    converge[i,s]<- sim$converge.success
    mean.parms$beta.a.ex=sim$out$mean$beta.a.ex
    mean.parms$beta.p.ex=sim$out$mean$beta.p.ex
    mean.parms$alpha.a=sim$out$mean$alpha.a
    mean.parms$alpha.p=sim$out$mean$alpha.p
    mean.parms$alpha.f=sim$out$mean$alpha.f
    sd.parms.init$vcovAban=sim$vcovAban
    sd.parms.init$vcovPred=sim$vcovPred
    sd.parms.init$alpha.f=sim$out$sd$alpha.f
    sd.parms.init$alpha.a=sim$out$sd$alpha.a
    sd.parms.init$alpha.p=sim$out$sd$alpha.p
    lambda.ex[i,,s] <- lambda.calc(parms=mean.parms, sd.parms=sd.parms.init, n.ex=1)$lambda
    lambda.un[i,,s] <- lambda.calc(parms=mean.parms, sd.parms=sd.parms.init, n.ex=0)$lambda
    cat("Finished",i,"of",n,"runs in ", s, "of", n.samples, "simulations") 
  }
}



#save####
save.image("~/Nest Survival Simulation/simulations/Tool_Sims.RData")
