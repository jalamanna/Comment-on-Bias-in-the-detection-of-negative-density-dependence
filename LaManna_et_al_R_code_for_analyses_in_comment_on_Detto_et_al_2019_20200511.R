### Analysis of CNDD in BCI and SERC data
### Code modified from Detto et al. 2019, data from Detto et al. 2019
### Code reproduces Fig. 1 in LaManna et al. comment on Detto et al. 2019
### by: J. LaManna, updated: 5-11-2020



#######################
### Load data
### First, download data files 'bci.mat' and 'serc.mat' from: https://github.com/mdetto/Bias-in-the-detection-of-negative-density-dependence
### Ensure that these two data files are placed into your R console's working directory prior to running the code below

library(R.matlab)

bci <- readMat("bci.mat")
bci = data.frame(dbh = bci$dbh, gx = bci$gx, gy = bci$gy, sp = unlist(bci$sp), treeID = bci$treeID)
serc <- readMat("serc.mat")
serc = data.frame(dbh = serc$dbh, gx = serc$gx, gy = serc$gy, sp = unlist(serc$sp), treeID = serc$treeID)


#######################
### Load functions

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

# 2.5th and 97.5th percentiles of mean dispersal distances simulated for a hypothetical 30 m max height tree species 
# These 100 replicates highlight the wide range of mean dispersal distances simulated for any given species in the Thomas-process null model
quantile(replicate(n = 100, mean.disp.dist(30)),c(0.025,0.975))

# DBH-height allometry function from ForestGEO CTFS Package (based on Chave 2005 Oecologia)
predht.asym=function(dbh,param)
{
 if(is.null(dim(param)))
  {
   ymax=param[1]
   a=param[2]
   b=param[3]
  }
 else
  {
   ymax=param[,1]
   a=param[,2]
   b=param[,3]
  }

 return(ymax*(1-exp(-a*dbh^b)))
}
htparam=c(41.7, 0.057, 0.748) # default CTFS height parameters from BCI data


DistWeighted = function(sp, dbh, gx, gy, Lx, Ly, DX, tr, quant, L, adult.mort, type) {
# quadrat grid
x = seq(0, Lx, by = DX)
y = seq(0, Ly, by = DX)
n = length(x) - 1
m = length(y) - 1

# inizialize variables
species = unique(sp)
S = length(species)
saplings = matrix(0, n*m, S)
adults = matrix(0, n*m, S)
adultDistWeighted = matrix(0, n*m, S)
BasalArea = rep(NA, times = S)
N = rep(NA, times = S)
Dcutoff = rep(NA, times = S)
sapquads = rep(NA, times = S)
adultquads = rep(NA, times = S)


for(s in 1:S) {
    use = which(sp == species[s] & dbh > 0)
    N[s] = length(use)
    BasalArea[s] = pi / 4 * sum((dbh[use] / 100) ^ 2) / (Lx * Ly / 10000)

    if(N[s] > tr) {
        x0 = gx[use]; y0 = gy[use]
        D = dbh[use]
        if(type == 'random') {D = D[sample(length(D))]}
        D50 = quantile(D, quant)	# select adult/sapling cutoff point
	  Dcutoff[s] = D50		# record DBH cutoff point
        
        # get adults and saplings position
        use1 = which(D <= D50)		# saplings
        use2 = which(D > D50)			# adults
        x1 = x0[use1]; y1 = y0[use1]	# sap
        x2 = x0[use2]; y2 = y0[use2]	# adults
	
	  # Wiegand-Moloney Thomas Process null model (fixes adult locations to preserve habitat heterogeneity/preferences, then disperses young away from adults given Clark's 2dT dispersal kernel)
		if(type == 'constant.disp.null') {
			alpha = clark.alpha(30)
			num.live = ceiling(length(use2)*(1-adult.mort)); num.dead = length(use2) - num.live
			repro.adults <- sample(1:length(use2), size = (length(use1) + num.dead), replace = T)
			saplocs = do.call('rbind', apply(data.frame(x2, y2)[repro.adults,], 1, dispersal.fun2, xlim = Lx, ylim = Ly, alpha = alpha))
			live.adults = sample(1:length(use2), size = num.live, replace = F)
			new.adults = sample(1:nrow(saplocs), size = num.dead, replace = F)
			x2 = c(x2[live.adults], saplocs[new.adults,1]); y2 = c(y2[live.adults], saplocs[new.adults,2])
			x1 = saplocs[which(1:nrow(saplocs) %in% new.adults == F),1]; y1 = saplocs[which(1:nrow(saplocs) %in% new.adults == F),2]						
		}
		if(type == 'allometric.disp.null') {
			max.height = predht.asym(dbh = max(D), param = htparam)
			mean.disp.distance = mean.disp.dist(max.height)	
			alpha = clark.alpha(mean.disp.distance)
			num.live = ceiling(length(use2)*(1-adult.mort)); num.dead = length(use2) - num.live
			repro.adults <- sample(1:length(use2), size = (length(use1) + num.dead), replace = T)
			saplocs = do.call('rbind', apply(data.frame(x2, y2)[repro.adults,], 1, dispersal.fun2, xlim = Lx, ylim = Ly, alpha = alpha))
			live.adults = sample(1:length(use2), size = num.live, replace = F)
			new.adults = sample(1:nrow(saplocs), size = num.dead, replace = F)
			x2 = c(x2[live.adults], saplocs[new.adults,1]); y2 = c(y2[live.adults], saplocs[new.adults,2])
			x1 = saplocs[which(1:nrow(saplocs) %in% new.adults == F),1]; y1 = saplocs[which(1:nrow(saplocs) %in% new.adults == F),2]
		}
        
        for(i in 1:n) {
            use1 = which(x1 >= x[i] & x1 < x[i] + DX)
            use2 = which(x2 >= x[i] & x2 < x[i] + DX)
            dx = x[i] - x2 + DX / 2
            for(j in 1:m) {
    		    LinInd = (i - 1) * m + j
                saplings[LinInd,s] = sum(y1[use1] >= y[j] & y1[use1] < y[j] + DX)
                adults[LinInd,s] = sum(y2[use2] >= y[j] & y2[use2] < y[j] + DX)
                dy = y[j] - y2 + DX / 2
                r2 = (dx * dx) + (dy * dy)
                adultDistWeighted[LinInd,s] = sum(1 / (pi * L ^ 2) * (L ^ 2 / (r2 + L ^ 2)) ^ 2) * DX ^ 2       
            }
        } 
	  sapquads[s] = sum(saplings[,s] > 0) 
	  adultquads[s] = sum(adults[,s] > 0) 
    }
}
result = list(saplings, adultDistWeighted, adults, species, N, BasalArea, Dcutoff, sapquads, adultquads)
names(result) = c("saplings", "adultDistWeighted", "adults", "species", "N", "BasalArea", "Dcutoff", "sapquads", "adultquads")
return(result)
}







###################################################################################################################
### RUN ANALYSIS
# parameters
quant = 0.5
L = 2 * 20 / pi; 		# mean dispersal = pi/2*L
DX = 20			# quadrat size
adult.mort = 0.10		# adult mortality (proportion of adults that die and are replaced by saplings in the dispersal-kernel model)
model = 'Ricker'		# model to estimate CNDD
null.type = 'allometric.disp.null'
null.its.k = 100		# iterations for the null model (final iterations will be this number times null.k)
null.k = 10

# Prepare BCI data
sp = bci$sp
bci$dbh = bci$dbh * 0.1
dbh = bci$dbh
gx = bci$gx
gy = bci$gy
species = unique(bci$sp)
S = length(species)
Lx = 1000; Ly = 500; bci.Area = 50
tr=25



## BCI analysis
set.seed(892)
bci.sl = rep(NA, times = S)
bci.slr = rep(NA, times = S)

bci2 = DistWeighted(sp, dbh, gx, gy, Lx, Ly, DX, tr, quant, L, adult.mort, type = 'observed')
bciran = DistWeighted(sp, dbh, gx, gy, Lx, Ly, DX, tr, quant, L, adult.mort, type = null.type)

# Observed and Thomas-process null model
for(i in 1:S) {
  if(bci2$N[i] > tr) {
    if(model == 'Ricker') {
      tbl = data.frame(A = bci2$adultDistWeighted[,i], S = bci2$saplings[,i], AR = bciran$adultDistWeighted[,i], SR = bciran$saplings[,i])
      m1 = glm(S ~ A + offset(log(A)), family = 'poisson', data = tbl)
      bci.sl[i] = coef(m1)[2]
      m2 = glm(SR ~ AR + offset(log(AR)), family = 'poisson', data = tbl)
      bci.slr[i] = coef(m2)[2]   
    } 
  }
}


# For repititions of the Thomas-process null model
bci.slr.mat = list()
for(z in 1:null.k) {
  rm(bciran)
  bci.slr.mat1 = matrix(NA, S, null.its.k)
  bciran = replicate(null.its.k, DistWeighted(sp, dbh, gx, gy, Lx, Ly, DX, tr, quant, L, adult.mort, type = null.type), simplify = F)
  for(j in 1:length(bciran)) {
    for(i in 1:S) {
      if(bciran[[j]]$N[i] > tr) {
        if(model == 'Ricker') {
          tbl = data.frame(AR = bciran[[j]]$adultDistWeighted[,i], SR = bciran[[j]]$saplings[,i])
          m2 = glm(SR ~ AR + offset(log(AR)), family = 'poisson', data = tbl)
          bci.slr.mat1[i,j] = coef(m2)[2]    
        }
      }
    }
  }
bci.slr.mat[[z]] = bci.slr.mat1
}
bci.slr.mat = do.call('cbind', bci.slr.mat)

bci.tpm = bci.slr.mat



# Relabel null model
null.type = 'random'
set.seed(892)
bci.slr = rep(NA, times = S)

# For repititions of the Relabel null model
bci.slr.mat = list()
for(z in 1:null.k) {
  rm(bciran)
  bci.slr.mat1 = matrix(NA, S, null.its.k)
  bciran = replicate(null.its.k, DistWeighted(sp, dbh, gx, gy, Lx, Ly, DX, tr, quant, L, adult.mort, type = null.type), simplify = F)
  for(j in 1:length(bciran)) {
    for(i in 1:S) {
      if(bciran[[j]]$N[i] > tr) {
        if(model == 'Ricker') {
          tbl = data.frame(AR = bciran[[j]]$adultDistWeighted[,i], SR = bciran[[j]]$saplings[,i])
          m2 = glm(SR ~ AR + offset(log(AR)), family = 'poisson', data = tbl)
          bci.slr.mat1[i,j] = coef(m2)[2]    
        }
      }
    }
  }
bci.slr.mat[[z]] = bci.slr.mat1
}
bci.slr.mat = do.call('cbind', bci.slr.mat)

bci.rlm = bci.slr.mat








####################################
# Prepare SERC data
sp = serc$sp
serc$dbh = serc$dbh * 0.1
dbh = serc$dbh
gx = serc$gx
gy = serc$gy
species = unique(serc$sp)
S = length(species)
Lx = 400; Ly = 400; serc.Area = 16
tr=11
null.type = 'allometric.disp.null'

## SERC analysis
set.seed(842)
serc.sl = rep(NA, times = S)
serc.slr = rep(NA, times = S)

serc2 = DistWeighted(sp, dbh, gx, gy, Lx, Ly, DX, tr, quant, L, adult.mort, type = 'observed')
sercran = DistWeighted(sp, dbh, gx, gy, Lx, Ly, DX, tr, quant, L, adult.mort, type = null.type)


# Observed and Thomas-process null model
for(i in 1:S) {
  if(serc2$N[i] > tr) {
    if(model == 'Ricker') {
      tbl = data.frame(A = serc2$adultDistWeighted[,i], S = serc2$saplings[,i], AR = sercran$adultDistWeighted[,i], SR = sercran$saplings[,i])
      m1 = glm(S ~ A + offset(log(A)), family = 'poisson', data = tbl)
      serc.sl[i] = coef(m1)[2]
      m2 = glm(SR ~ AR + offset(log(AR)), family = 'poisson', data = tbl)
      serc.slr[i] = coef(m2)[2]   
    } 
  }
}


# For repititions of the Thomas-process null model
serc.slr.mat = list()
for(z in 1:null.k) {
  rm(sercran)
  serc.slr.mat1 = matrix(NA, S, null.its.k)
  sercran = replicate(null.its.k, DistWeighted(sp, dbh, gx, gy, Lx, Ly, DX, tr, quant, L, adult.mort, type = null.type), simplify = F)
  for(j in 1:length(sercran)) {
    for(i in 1:S) {
      if(sercran[[j]]$N[i] > tr) {
        if(model == 'Ricker') {
          tbl = data.frame(AR = sercran[[j]]$adultDistWeighted[,i], SR = sercran[[j]]$saplings[,i])
          m2 = glm(SR ~ AR + offset(log(AR)), family = 'poisson', data = tbl)
          serc.slr.mat1[i,j] = coef(m2)[2]   
        } 
      }
    }
  }
serc.slr.mat[[z]] = serc.slr.mat1
}
serc.slr.mat = do.call('cbind', serc.slr.mat)

serc.tpm = serc.slr.mat



# Relabel null model
null.type = 'random'
set.seed(842)
serc.slr = rep(NA, times = S)

# For repititions of the Relabel null model
serc.slr.mat = list()
for(z in 1:null.k) {
  rm(sercran)
  serc.slr.mat1 = matrix(NA, S, null.its.k)
  sercran = replicate(null.its.k, DistWeighted(sp, dbh, gx, gy, Lx, Ly, DX, tr, quant, L, adult.mort, type = null.type), simplify = F)
  for(j in 1:length(sercran)) {
    for(i in 1:S) {
      if(sercran[[j]]$N[i] > tr) {
        if(model == 'Ricker') {
          tbl = data.frame(AR = sercran[[j]]$adultDistWeighted[,i], SR = sercran[[j]]$saplings[,i])
          m2 = glm(SR ~ AR + offset(log(AR)), family = 'poisson', data = tbl)
          serc.slr.mat1[i,j] = coef(m2)[2]   
        } 
      }
    }
  }
serc.slr.mat[[z]] = serc.slr.mat1
}
serc.slr.mat = do.call('cbind', serc.slr.mat)

serc.rlm = serc.slr.mat







##################################
### PLOT RESULTS



bci.random.label = bci.rlm
serc.random.label = serc.rlm
bci.allometric = bci.tpm
serc.allometric = serc.tpm
bci.sl2 = bci.sl[!is.na(bci.sl)]
serc.sl2 = serc.sl[!is.na(serc.sl)]
bci.ba = bci2$BasalArea[!is.na(bci.sl)]
serc.ba = serc2$BasalArea[!is.na(serc.sl)]





### Figures 1A and 1B
par(mfrow = c(2,1))

# BCI All species
bci.null.rl.medCNDD = apply(bci.random.label, 2, median, na.rm = T)
(bci.obs.median = median(bci.sl, na.rm = T))
serc.null.rl.medCNDD = apply(serc.random.label, 2, median, na.rm = T)
(serc.obs.median = median(serc.sl, na.rm = T))
bci.null.allo.medCNDD = apply(bci.allometric, 2, median, na.rm = T)
serc.null.allo.medCNDD = apply(serc.allometric, 2, median, na.rm = T)
hist(bci.null.rl.medCNDD, xlim = c(-2.5, 0), breaks = 20, las = 1, col = "gray85",
	xlab = "Median CNDD at BCI (all species)", main = "")
hist(bci.null.allo.medCNDD, breaks = 30, las = 1, col = "lightblue", add = T)
abline(v = bci.obs.median, col = "red", lwd = 3); abline(v = 0, lty = 2)
sum(bci.null.rl.medCNDD < bci.obs.median) / length(bci.null.rl.medCNDD)
sum(bci.null.allo.medCNDD < bci.obs.median) / length(bci.null.allo.medCNDD)

# SERC All species
hist(serc.null.rl.medCNDD, xlim = c(-2.5, 0.0), breaks = 5, las = 1, col = "gray85",
	xlab = "Median CNDD at SERC (all species)", main = "")
hist(serc.null.allo.medCNDD, breaks = 30, las = 1, col = "lightblue", add = T)
abline(v = serc.obs.median, col = "red", lwd = 3); abline(v = 0, lty = 2)
sum(serc.null.rl.medCNDD < serc.obs.median) / length(serc.null.rl.medCNDD)
sum(serc.null.allo.medCNDD > serc.obs.median) / length(serc.null.allo.medCNDD)




### Figures 1C and 1D
par(mfrow = c(2,1))

# All species
bci.null.rl.medCNDD = apply(bci.random.label, 2, median, na.rm = T)
(bci.obs.median = median(bci.sl, na.rm = T))
serc.null.rl.medCNDD = apply(serc.random.label, 2, median, na.rm = T)
(serc.obs.median = median(serc.sl, na.rm = T))
bci.null.allo.medCNDD = apply(bci.allometric, 2, median, na.rm = T)
serc.null.allo.medCNDD = apply(serc.allometric, 2, median, na.rm = T)
(obs.diff = median(bci.sl, na.rm = T) - median(serc.sl, na.rm = T))
null.diff.rl = bci.null.rl.medCNDD - serc.null.rl.medCNDD
null.diff.allo = bci.null.allo.medCNDD - serc.null.allo.medCNDD
hist(null.diff.rl, xlim = c(-2.5, 0), breaks = 20, las = 1, col = "gray85",
	xlab = "Difference in median CNDD (all species)", main = "")
hist(null.diff.allo, breaks = 50, las = 1, col = "lightblue", add = T)
abline(v = obs.diff, col = "red", lwd = 3); abline(v = 0, lty = 2)
sum(null.diff.rl < obs.diff) / length(null.diff.rl)

# Rare species
bci.null.rl.medCNDD = apply(bci.random.label[which(bci2$BasalArea <= 0.1),], 2, median, na.rm = T)
(bci.obs.median = median(bci.sl[which(bci2$BasalArea <= 0.1)], na.rm = T))
serc.null.rl.medCNDD = apply(serc.random.label[which(serc2$BasalArea <= 0.1),], 2, median, na.rm = T)
(serc.obs.median = median(serc.sl[which(serc2$BasalArea <= 0.1)], na.rm = T))
bci.null.allo.medCNDD = apply(bci.allometric[which(bci2$BasalArea <= 0.1),], 2, median, na.rm = T)
serc.null.allo.medCNDD = apply(serc.allometric[which(serc2$BasalArea <= 0.1),], 2, median, na.rm = T)
(obs.diff = median(bci.sl[which(bci2$BasalArea <= 0.1)], na.rm = T) - median(serc.sl[which(serc2$BasalArea <= 0.1)], na.rm = T))
null.diff.rl = bci.null.rl.medCNDD - serc.null.rl.medCNDD
null.diff.allo = bci.null.allo.medCNDD - serc.null.allo.medCNDD
hist(null.diff.rl, xlim = c(-2.5, 0.0), breaks = 30, las = 1, col = "gray85",
	xlab = "Difference in median CNDD (rare species)", main = "")
hist(null.diff.allo, breaks = 50,  las = 1, col = "lightblue", add = T)
abline(v = obs.diff, col = "red", lwd = 3); abline(v = 0, lty = 2)
sum(null.diff.rl < obs.diff) / length(null.diff.rl)

































