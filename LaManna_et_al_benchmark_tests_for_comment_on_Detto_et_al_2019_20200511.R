### Benchmark tests presented in Fig. 2 of LaManna et al. comment on Detto et al. 2019 
### showing the random-label model used by Detto et al. induces false density dependence
###
### R script to simulate data with known value of CNDD, simulated error,
### and dispersal of seeds outside the forest plot for CTFS-ForestGEO analyses.
###
### By: J. A. LaManna, updated 04-16-2020


# Load required packages
library(doBy)
library(reshape)
library(VGAM)
library(emdbook)
library(nlme)
library(ggplot2)
library(RColorBrewer)
library(mvtnorm)
library(doBy)
library(ggplot2)
library(pscl)
library(MASS)
library(boot)
library(spatstat)
library(vegan)



################################################################################################
# I. Functions
################################################################################################
### Functions to simulate data with known CNDD
# n = number of quadrats in simulated forest plot
# meanTrees = number of mean adult trees per quadrat
# lambda = per capita recruitment rate (in absence of density dependence)
# trueCNDD = conspecific density dependence (0 = no CNDD)
# theta = negative binomial overdispersion parameter
# d = proportion of seeds dispersing outside of the simulated forest plot

################################################################################################
### Fit Functions

### Load functions for analyses

fit.ricker.cndd.dist.wgt = function(data){
AR = data$distwgtTrees; SR = data$recruits
m2 = glm(SR ~ AR + offset(log(AR)), family = 'poisson')
return(m2)
}


#############################################################################################################
### Other necessary functions

seed = function(r, alpha=100) {(1/(pi*alpha*((1+((r^2)/alpha))^2)))}	# Clark's 2Dt p = 1 

# Function to calculate distance-weighted adult abundances for each quadrat
dist.weighted.abund.function = function(test, quads) {
adultlocs = test[,c("gx", "gy")]
quadlocs = quads[,c("cx", "cy")]+10
xdist = crossdist(adultlocs[,1], adultlocs[,2], quadlocs[,1], quadlocs[,2])
xdist.weighted = seed(xdist, alpha = exp(6)) * (1/seed(0.1, alpha = exp(6)))
return(apply(xdist.weighted, 2, sum))
}

# Clark's 2Dt seed dispersal kernel
dispersal.fun = function(x, y, alpha, n) {
seedp = function(r, alpha=100) {(1/(pi*alpha*((1+((r^2)/alpha))^2)))*2*pi*r}	# Clark's 2Dt p = 1 (integrated in 2 polar dimensions)
disp.dists = seq(0, 10000, by = 1)
rho <- sample(disp.dists, size = n, replace = TRUE, prob = seedp(r = disp.dists, alpha = alpha))
theta <- runif(n, 0, 2*pi)
x1 <- (rho * cos(theta)) + x
y1 <- (rho * sin(theta)) + y
result = data.frame(x1, y1)
return(result)
}

dispersal.fun2 = function(loc, xlim = plotwidth, ylim = plotheight, alpha) {	# Recruits disperse across plot edges in a torus
test = dispersal.fun(loc[1], loc[2], alpha = alpha, n = 1)
test[,1] = test[,1] %% xlim		# Torus
test[,2] = test[,2] %% ylim
return(data.frame(x1 = test[,1], y1 = test[,2]))
}

# Thomson et al 2011 J. of Ecol. Relationship between max tree height and mean dispersal distance (plus observed error around regression fit)
mean.disp.dist = function(max.height) {10 ^ ((0.1975 + (0.936 * log10(max.height))) + rnorm(length(max.height), mean = 0, sd = 0.245))}		# Relationship between maximum tree height and mean dispersal distance with observed error from Thomson et al 2011 J. of Ecol.
clark.alpha = function(mean.disp, shape = 1) {exp(log(mean.disp) - (0.5*log(pi) + log(gamma(shape - 1/2)) - log(2) - log(gamma(shape))))^2}	# Calculate alpha from Clark's 2Dt model given mean dispersal distance

# Function to disperse saplings around fixed adult locations given a mean dispersal distance (Clark's 2dT kernel)
dispersal.kernel.function = function(mean.disp = 30, adultlist, saplinglist, Lx, Ly) {
	x2 = adultlist$gx; y2 = adultlist$gy
	alpha = clark.alpha(mean.disp)
	repro.adults <- sample(1:nrow(adultlist), size = nrow(saplinglist), replace = T)
	saplocs = do.call('rbind', apply(data.frame(x2, y2)[repro.adults,], 1, dispersal.fun2, xlim = Lx, ylim = Ly, alpha = alpha))
	return(saplocs)					
}

# Function to ID quadrat given a location
quad.ID.function = function(Lx, Ly, DX, x1, y1, S = 1) {
x = seq(0, Lx, by = DX)
y = seq(0, Ly, by = DX)
n = length(x) - 1
m = length(y) - 1
saplings = matrix(0, n*m, S)
    for(i in 1:n) {
        use1 = which(x1 >= x[i] & x1 < x[i] + DX)
        for(j in 1:m) {
	    LinInd = (i - 1) * m + j
          saplings[LinInd] = sum(y1[use1] >= y[j] & y1[use1] < y[j] + DX)
        }
    } 
return(saplings)
}


################################################################################################
### Data simulation functions

# Simple Model (no error)
sim.data.simple.ricker = function(h = 25, w = 50, meanTrees, r, trueCNDD, theta){
n = h * w
trueTrees = rnbinom(n, size = theta, mu = meanTrees)
cy1 = (1:h * 20) - 20
cx1 = (1:w * 20) - 20
quad.centers = expand.grid(cy1, cx1); quad.centers = data.frame(cx = quad.centers[,2], cy = quad.centers[,1])
data1 = data.frame(trueTrees, cx = quad.centers$cx, cy = quad.centers$cy)
data1$quadrat = c(1:nrow(data1))
trees.expanded = data1[rep(seq.int(1,nrow(data1)), data1$trueTrees),]
trees.expanded$gx = trees.expanded$cx + runif(nrow(trees.expanded),0.001,19.999)
trees.expanded$gy = trees.expanded$cy + runif(nrow(trees.expanded),0.001,19.999)
data1$dist.wgt.adults = dist.weighted.abund.function(trees.expanded, data1)
recruits = data1$dist.wgt.adults * exp(r + (trueCNDD * data1$dist.wgt.adults)) 	# Ricker population model
data=data.frame(consppTrees = trueTrees, distwgtTrees = data1$dist.wgt.adults, recruits)
}


# Error Model
# same thing but add process error to observed recruits per quadrat, simulates demographic stochasticity
sim.data.error.ricker = function(h = 25, w = 50, meanTrees, r, trueCNDD, theta){
n = h * w
trueTrees = rnbinom(n, size = theta, mu = meanTrees)
cy1 = (1:h * 20) - 20
cx1 = (1:w * 20) - 20
quad.centers = expand.grid(cy1, cx1); quad.centers = data.frame(cx = quad.centers[,2], cy = quad.centers[,1])
data1 = data.frame(trueTrees, cx = quad.centers$cx, cy = quad.centers$cy)
data1$quadrat = c(1:nrow(data1))
trees.expanded = data1[rep(seq.int(1,nrow(data1)), data1$trueTrees),]
trees.expanded$gx = trees.expanded$cx + runif(nrow(trees.expanded),0.001,19.999)
trees.expanded$gy = trees.expanded$cy + runif(nrow(trees.expanded),0.001,19.999)
data1$dist.wgt.adults = dist.weighted.abund.function(trees.expanded, data1)
truerecruits = data1$dist.wgt.adults * exp(r + (trueCNDD * data1$dist.wgt.adults)) 	# Ricker population model
recruits = rnbinom(n, size = theta, mu = truerecruits)
data=data.frame(consppTrees = trueTrees, truerecruits, distwgtTrees = data1$dist.wgt.adults, recruits)
data2 = data.frame(recruits, cx = quad.centers$cx, cy = quad.centers$cy)
data2$quadrat = c(1:nrow(data2))
saps.expanded = data2[rep(seq.int(1,nrow(data2)), data2$recruits),]
saps.expanded$gx = saps.expanded$cx + runif(nrow(saps.expanded),0.001,19.999)
saps.expanded$gy = saps.expanded$cy + runif(nrow(saps.expanded),0.001,19.999)
saps.expanded$sap = c("sap")
trees.expanded$sap = c("adult")
inds = rbind(trees.expanded[,-1], saps.expanded[,-1])
result = list(data, inds, data1)
names(result) = c("data", "inds", "quads"); return(result)
}


# Dispersal and error model
# same as error version but some fraction (d) of recruits globally dispersed
sim.data.dispersal.ricker = function(h = 25, w = 50, meanTrees, r, trueCNDD, theta, d){
n = h * w
trueTrees = rnbinom(n, size = theta, mu = meanTrees)
cy1 = (1:h * 20) - 20
cx1 = (1:w * 20) - 20
quad.centers = expand.grid(cy1, cx1); quad.centers = data.frame(cx = quad.centers[,2], cy = quad.centers[,1])
data1 = data.frame(trueTrees, cx = quad.centers$cx, cy = quad.centers$cy)
data1$quadrat = c(1:nrow(data1))
trees.expanded = data1[rep(seq.int(1,nrow(data1)), data1$trueTrees),]
trees.expanded$gx = trees.expanded$cx + runif(nrow(trees.expanded),0.001,19.999)
trees.expanded$gy = trees.expanded$cy + runif(nrow(trees.expanded),0.001,19.999)
data1$dist.wgt.adults = dist.weighted.abund.function(trees.expanded, data1)
totRecruits = data1$dist.wgt.adults * exp(r + (trueCNDD * data1$dist.wgt.adults)) 	# Ricker population model (with dist. weighted adult abundances)
localRecruits = totRecruits * (1 - d)
recruits = localRecruits + ((sum(totRecruits) * d) / n)
recruits = rnbinom(n, size = theta, mu = recruits)
data=data.frame(consppTrees = trueTrees, truerecruits = totRecruits, distwgtTrees = data1$dist.wgt.adults, recruits)
data2 = data.frame(recruits, cx = quad.centers$cx, cy = quad.centers$cy)
data2$quadrat = c(1:nrow(data2))
saps.expanded = data2[rep(seq.int(1,nrow(data2)), data2$recruits),]
saps.expanded$gx = saps.expanded$cx + runif(nrow(saps.expanded),0.001,19.999)
saps.expanded$gy = saps.expanded$cy + runif(nrow(saps.expanded),0.001,19.999)
saps.expanded$sap = c("sap")
trees.expanded$sap = c("adult")
inds = rbind(trees.expanded[,-1], saps.expanded[,-1])
result = list(data, inds, data1)
names(result) = c("data", "inds", "quads"); return(result)
}



# Dispersal and error model with adults mortality
# same as error/dispersal version but adds measurement error to adults
sim.data.error.dispersal.ricker = function(h = 25, w = 50, meanTrees, r, trueCNDD, theta, d, m){
n = h * w
trueTrees = rnbinom(n, size = theta, mu = meanTrees)
cy1 = (1:h * 20) - 20
cx1 = (1:w * 20) - 20
quad.centers = expand.grid(cy1, cx1); quad.centers = data.frame(cx = quad.centers[,2], cy = quad.centers[,1])
data1 = data.frame(trueTrees, cx = quad.centers$cx, cy = quad.centers$cy)
data1$quadrat = c(1:nrow(data1))
trees.expanded = data1[rep(seq.int(1,nrow(data1)), data1$trueTrees),]
trees.expanded$gx = trees.expanded$cx + runif(nrow(trees.expanded),0.001,19.999)
trees.expanded$gy = trees.expanded$cy + runif(nrow(trees.expanded),0.001,19.999)
data1$dist.wgt.adults = dist.weighted.abund.function(trees.expanded, data1)
totRecruits = data1$dist.wgt.adults * exp(r + (trueCNDD * data1$dist.wgt.adults)) 	# Ricker population model (with dist. weighted adult abundances)
localRecruits = totRecruits * (1 - d)
recruits = localRecruits + ((sum(totRecruits) * d) / n)
recruits = rnbinom(n, size = theta, mu = recruits)
TRUE.dist.weights = data1$dist.wgt.adults
num.live = ceiling(sum(trueTrees)*(1-m))
trees.expandedR = trees.expanded[sample(1:nrow(trees.expanded), size = num.live, replace = F),]
data2 = unique(trees.expandedR[,c("cx", "cy", "quadrat")])
OBS.dist.weights = dist.weighted.abund.function(trees.expandedR, data1)
data=data.frame(consppTrees = trueTrees, truerecruits = totRecruits, distwgtTrees = OBS.dist.weights, recruits)
data2 = data.frame(recruits, cx = quad.centers$cx, cy = quad.centers$cy)
data2$quadrat = c(1:nrow(data2))
saps.expanded = data2[rep(seq.int(1,nrow(data2)), data2$recruits),]
saps.expanded$gx = saps.expanded$cx + runif(nrow(saps.expanded),0.001,19.999)
saps.expanded$gy = saps.expanded$cy + runif(nrow(saps.expanded),0.001,19.999)
saps.expanded$sap = c("sap")
trees.expandedR$sap = c("adult")
inds = rbind(trees.expandedR[,-1], saps.expanded[,-1])
result = list(data, inds, data1)
names(result) = c("data", "inds", "quads"); return(result)
}




################################################################################################
# I. Simple: no observation or measurement error, no dispersal
################################################################################################
# Ricker models that do not use distance-weighted adult abundance are blind to number of 
# adults trees in neighboring quadrats

# Simple Ricker Simulation
h = 25
w = 50
meanTrees = 2
trueCNDD = -0.1	# True conspecific negative density-dependence (0 = no CNDD; lower values = stronger CNDD)
r = 0
theta = 1
data = sim.data.simple.ricker(h=h,w=w,meanTrees,r=r,trueCNDD=trueCNDD,theta=theta)
m = fit.ricker.cndd.dist.wgt(data) 
summary(m)

plot(data$distwgtTrees, data$recruits, ylim=c(0,10))
x = seq(min(data$distwgtTrees), max(data$distwgtTrees), length = 1000)
f1 = function(x) {x*exp(coef(m)[1] + (x*coef(m)[2]))}
lines(x, f1(x))




################################################################################################
# II. Add process error to recruits
################################################################################################
# add process error to number of recruits in each quadrat

meanTrees = 1
h = 25
w = 50
trueCNDD = -0.10	# True conspecific negative density-dependence (0 = no CNDD; lower values = stronger CNDD)
r = 0.00
theta = 10
data = sim.data.error.ricker(h = h, w = w, meanTrees = meanTrees, r = r, trueCNDD = trueCNDD, theta = theta)
m = fit.ricker.cndd.dist.wgt(data$data) 
summary(m)

plot(data$data$distwgtTrees, data$data$recruits, ylim=c(0,10))
x = seq(min(data$data$distwgtTrees), max(data$data$distwgtTrees), length = 1000)
f1 = function(x) {x*exp(coef(m)[1] + (x*coef(m)[2]))}
lines(x, f1(x))




################################################################################################
# III. Add dispersal and process error
################################################################################################
# Incorporate dispersal out of the plot, ((1-d)*100) recruits stay put, (d*100) are globally dispersed

meanTrees = 0.145	# Mean number of adult trees per quadrat
h = 25		# Number of quadrats the forest plot is in height
w = 50		# Number of quadrats the forest plot is in width
trueCNDD = 0	# True conspecific negative density-dependence (0 = no CNDD; lower values = stronger CNDD)
r = 0.5		# Density-independent population growth rate 
theta = 2		# Error
d = 0.05		# Dispersal Factor
data = sim.data.dispersal.ricker(h = h, w = w, meanTrees = meanTrees, r = r, trueCNDD = trueCNDD, theta = theta, d = d)
m = fit.ricker.cndd.dist.wgt(data$data) 
summary(m)

plot(data$data$distwgtTrees, data$data$recruits, ylim=c(0,10))
x = seq(min(data$data$distwgtTrees), max(data$data$distwgtTrees), length = 1000)
f1 = function(x) {x*exp(coef(m)[1] + (x*coef(m)[2]))}
lines(x, f1(x))



################################################################################################
# IV. Add dispersal, process error, and measurement error
################################################################################################
# Incorporate dispersal out of the plot, ((1-d)*100) recruits stay put, (d*100) are globally dispersed

set.seed(256)

meanTrees = 0.043	# Mean number of adult trees per quadrat
h = 25		# Number of quadrats the forest plot is in height
w = 50		# Number of quadrats the forest plot is in width
trueCNDD = 0	# True conspecific negative density-dependence in recruitment (0 = no CNDD; lower values = stronger CNDD)
r = -0.75		# Density-independent recruitment 
theta = 2		# Error
d = 0.05		# Dispersal Factor
m = 0.10		# Adult Mortality
its = 1000		# Number of iterations for null model
data = sim.data.error.dispersal.ricker(h = h, w = w, meanTrees = meanTrees, r = r, trueCNDD = trueCNDD, theta = theta, d = d, m = m)
m = fit.ricker.cndd.dist.wgt(data$data) 
summary(m)

plot(data$data$distwgtTrees, data$data$recruits, ylim=c(0,10))
x = seq(min(data$data$distwgtTrees), max(data$data$distwgtTrees), length = 1000)
f1 = function(x) {x*exp(coef(m)[1] + (x*coef(m)[2]))}
lines(x, f1(x))


# Random-label Null Model
rlm.r = c()
rlm.cndd = c()
for(i in 1:its) {
test = data$inds
quadsR = data$quads
test$sapR = test$sap[sample(1:nrow(test), replace = F)]
test.ad = test[which(test$sapR == "adult"),]
quadsR$dist.wgts = dist.weighted.abund.function(test.ad, data$quads)
test$value = c(1)
test2 = aggregate(value ~ quadrat, FUN = length, data = test[which(test$sapR == "sap"),])
test3 = merge(data$quads, test2, by = "quadrat", all.x = T)
test3$value[is.na(test3$value)] = 0
rlm.data = data.frame(distwgtTrees = quadsR$dist.wgts, recruits = test3$value)
mR = fit.ricker.cndd.dist.wgt(rlm.data) 
rlm.r[i] = coef(mR)[1]
rlm.cndd[i] = coef(mR)[2]
}

# Dispersal-kernel Null Model
dkm.r = c()
dkm.cndd = c()
for(i in 1:its) {
test = data$inds
quadsR = data$quads
test.ad = test[which(test$sap == "adult"),]
test.sap = test[which(test$sap == "sap"),]
saplocs = dispersal.kernel.function(mean.disp = 30, adultlist = test.ad, saplinglist = test.sap, Lx = w*20, Ly = h*20)
quadsR$value = quad.ID.function(Lx = w*20, Ly = h*20, DX = 20, x1 = saplocs[,1], y1 = saplocs[,2])
dkm.data = data.frame(distwgtTrees = quadsR$dist.wgt.adults, recruits = quadsR$value)
mR = fit.ricker.cndd.dist.wgt(dkm.data) 
dkm.r[i] = coef(mR)[1]
dkm.cndd[i] = coef(mR)[2]
}




# Plots showing the range of estimates produced by the null model, the truth (blue), and observed CNDD (red):
par(mfrow=c(2,1))
hist(rlm.cndd, xlim = range(c(rlm.cndd, dkm.cndd, coef(m)[2], trueCNDD))*1.05, ylim = c(0,70), breaks = 60, col = "gray85", main = "Random-label vs. Dispersal-kernel Null Model", xlab = "CNDD in recruitment", las = 1)
hist(dkm.cndd, xlim = range(c(rlm.cndd, dkm.cndd, coef(m)[2], trueCNDD))*1.05, breaks = 60, col = "lightblue", las = 1, add = T)
hist(rlm.cndd, xlim = range(c(rlm.cndd, dkm.cndd, coef(m)[2], trueCNDD))*1.05, breaks = 60, col = "gray85", las = 1, add = T)
hist(dkm.cndd, xlim = range(c(rlm.cndd, dkm.cndd, coef(m)[2], trueCNDD))*1.05, breaks = 60, col = rgb(t(col2rgb(c("lightblue"))), alpha=100, maxColorValue = 255), las = 1, add = T)
abline(v = trueCNDD, col = "blue", lwd = 3); abline(v = coef(m)[2], col = "red", lwd = 3)
















