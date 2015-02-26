#Code for Hamm and Fordyce 201X - Evolution

#This file contains all the data and code necessary to reproduce the analyses from our paper. I've commented out the lines for reading in the data because all the objects are saved in the R file. 

#This code was run on a 15" MBP with a 2.6GHz Intel Core i7 processor that had 16GB of RAM running OSX 10.9.4, and R 3.1.1.

set.seed(8)
setwd("") # set your working directory
# save(list=ls(), file = "Hamm_Fordyce_data_1.R")
load("Hamm_Fordyce_data_1.R") # here is an R file that contains all of the objects below in case you don't want to run all the code

library(phytools) # v0.3-91

# read in the tree
# Nymphalidae <- read.tree("Wahl_Beast_New4.tre")

# read in the covariates
# Covs <- read.csv("Nym_covariates_2.csv", row.names = 1)

# prune the tree and data such that only taxa present in both are included in analyses
Nym.pruned <- treedata(phy = Nymphalidae, data = Covs, sort = TRUE, warnings = TRUE)
ngenera <- length(Nym.pruned$phy$tip.label)

# Sort names to match, use drop=F to prevent class conversion when only one column is present
Covs <- Covs[Nym.pruned$phy$tip.label, , drop = FALSE]

# make objects with variable of interest
Hosts <- Covs[, 1] # The host breadth, which is the number of plant families a genus feeds upon
names(Hosts) <- row.names(Covs)
Hosts

Nsp <- Covs[, 2] # The number of species in a genus of Nymphalidae
names(Nsp) <- row.names(Covs)
Nsp

###################################################
# detecting rate heterogeneity in the tree
# using iteRates for the parametric rates comparison
library(iteRates) #  v3.1
trimmed_tree <- trimTree(Nym.pruned$phy, Time = 20) # trimming the last 20 mya off the tree to meet the assumptions of the PRC. Below are the function calls for: exponential, Weibull, rate variable, and a call for all models. 

Nym.exp <- comp.subs(trimmed_tree$t.tree, mod.id = c(1, 0, 0, 0))

Nym.Wei <- comp.subs(trimmed_tree$t.tree, mod.id = c(0, 1, 0, 0))

Nym.rv <- comp.subs(trimmed_tree$t.tree, mod.id = c(0, 0, 0, 1))

Nym.mod <- comp.subs(trimmed_tree$t.tree, mod.id = c(1, 1, 0, 1))

par(mfrow=c(1, 3)) # plot the trees to see what nodes stick out
color.tree.plot(Nym.exp, trimmed_tree$t.tree, p.thres = 0.05, cex = 0.3, label.offset = 1, use.edge.length = TRUE, edge.width = 1.5, font = 3, no.margin = TRUE, root.edge = TRUE)

color.tree.plot(Nym.Wei, trimmed_tree$t.tree, p.thres = 0.05, cex = 0.3, label.offset = 1, use.edge.length = TRUE, edge.width = 1.5, font = 3, no.margin = TRUE, root.edge = TRUE)

color.tree.plot(Nym.rv, trimmed_tree$t.tree, p.thres = 0.05, cex = 0.3, label.offset = 1, use.edge.length = TRUE, edge.width = 1.5, font = 3, no.margin = TRUE, root.edge = TRUE)

par(mfrow = c( 1, 1))
# pdf(file="Nym_iteRates.pdf", bg="white")
color.tree.plot(Nym.mod, trimmed_tree$t.tree, p.thres = 0.05, cex = 0.3, label.offset = 1, use.edge.length = TRUE, edge.width = 1.5, font = 3, no.margin = TRUE, root.edge = TRUE)
# dev.off()

# get oldest node
tail(sort(branching.times(Nymphalidae))) # oldest node if 104 mya. 


# turboMEDUSA here was run in the terminal to take advantage of parallel processiing, you may not be able to run the parallel version of this code in the R GUI.
n.taxa <- as.vector(as.numeric(Nym.pruned$data[, 2]))
taxon <- rownames(Nym.pruned$data)
Mrich <- as.data.frame(cbind(taxon, n.taxa))
Mrich
library(turboMEDUSA) # v0.12
tMED1 <- runTurboMEDUSA(phy = Nym.pruned$phy, richness = Mrich, model.limit = 20, stop = "model.limit", criterion = "aicc", mc = TRUE, num.cores = 6, shiftCut = "both")
summarize.turboMEDUSA(tMED1, time = FALSE, cex = 0.2, plotTree = FALSE)

# using the model.limit or threshold gives the same answer for breakpoints
tMED2 <- runTurboMEDUSA(phy = Nym.pruned$phy, richness = Mrich, model.limit = 20, stop = "threshold", criterion = "aicc", mc = TRUE, num.cores = 6)
summarize.turboMEDUSA(tMED2, time = FALSE, cex = 0.2, plotTree = FALSE)


# Taking the inferred net diversification and relative extinction rates and converting them into the per-lineage speciation and extinction rates
getBD <- function (r, epsilon = NA) {
 	if (is.na(epsilon)) {epsilon <- 0} 
 	b <- r / (1-epsilon)
 	d <- b-r  
 	return(list( b= b, d = d))
}

# split node 1
BD1 <- getBD(0.059212, 9.0977e-01)
# split node 2
BD2 <- getBD(0.184138, 7.5193e-01)
# split node 3
BD3 <- getBD(0.100298, 6.7908e-01)
# split node 4
BD4 <- getBD(0.029644, 5.6883e-08)
# split node 5
BD5 <- getBD(0.030351, 2.2954e-07)
BD_rates_tMED1 <- matrix(data = NA, 4, 2, dimnames = list(c("Split node 2", "Split node 3", "Split node 4", "Split node 5"), c("lambda", "mu")))
BD_rates_tMED1[, 1] <- c(BD2$b, BD3$b, BD4$b, BD5$b)
BD_rates_tMED1[, 2] <- c(BD2$d, BD3$d, BD4$d, BD5$d)
BD_rates_tMED1 # note in the paper we do not mark node 1, so everything is shifted forward one node. 

# can eyeball and see the daughters of node 2 had an increase in diversification rate. For node 3 it's less clear. The estimated BAMM rates are below. 
#Node 2:
# Dircennina & Godyridina 
0.74228242 - 5.581444e-01
# Ithomiini (sister)
0.65623407 - 5.970221e-01

# Node 3:
# Satyrinae diversification rate
0.31253272 - 2.122347e-01
# Charaxine (sister)
0.65623407 - 5.970221e-01

# Node 4
# Kalimini
0.029644 - 5.6883e-08
# Junoniini (sister)
0.65623407 - 5.970221e-01

# Node 5 part of Satyrini
0.03035 - 6.967e-09
# other part (sister)
0.31253272 - 2.122347e-01

# BAMM analysis - nb that the BAMM program was run the terminal and analyzed here. Disclaimer, I followed a lot of the excellent BAMM website for code and tricks http://bamm-project.org/. You should check it out too!

library(BAMMtools) # v2.0.2
library(coda) # v0.16.1

# here we generate some priors for our data set, it outputs a text file to the working directory
# setBAMMpriors(Nym.pruned$phy, total.taxa = 500)

# read in the mcmc file
mcmc1 <- read.table("mcmc_out.txt", sep=",", header=TRUE)

# look at generations to be pruned for burnin
plot(mcmc1[1:1251, 1], mcmc1[1:1251, 4], pch = 19, las = 1) # looks like it has reached stationarity after ~10%

#First thing we will do is look at the autocorrelation and effective sample size in our data
# autocorr(as.mcmc(mcmc1), lags=c(0, 1, 5, 10, 50), relative = TRUE) # going to use the lag = 50
autocorr.plot(mcmc1, lag.max = 50) # autocorrelation seems to drop off pretty quickly for N_shifts, event_rate and accept_rate. 

effectiveSize(mcmc1) # calculate the effective sample size using coda

burnstart <- floor(0.1 * nrow(mcmc1))
postburn <- mcmc1[burnstart:nrow(mcmc1), ] # remove the burnin
plot(postburn[1:1127, 1], postburn[1:1127, 4], pch = 19, las = 1) # looks good, no pattern!

#check autocorrelation of postburn data
# autocorr(as.mcmc(postburn), lags=c(0, 1, 5, 10, 50), relative = TRUE) # note how rapidly autocorrelation goes away once the burnin is removed
effectiveSize(postburn) # Note that by removing the burin we significantly increased the effective size of our samples
autocorr.plot(postburn, lag.max = 50) # with the brunin removed the autocorrelation is effectively gone

# How many rate shifts? Look at postior odds ratios
post.probs <- table(postburn[1:1127, 2]) / nrow(postburn)
names(post.probs)
post.probs["2"] / post.probs["3"] # 2.23 posterior odds ratio, this means that 2 shifts in the tree are 2.24 times more likely to have occured than 3 shifts
post.probs["2"] / post.probs["4"] # 8.25
post.probs["2"] / post.probs["5"] # 64.54


# read in the event data from the BAMM run
event_data <- read.table("event_data.txt", sep = ",", header = TRUE)
ev1 <- getEventData(Nym.pruned$phy, eventdata = event_data, burnin = 0.1, type = "diversification")
Nym.shift.probs <- summary(ev1) # 63% of postior had two shifts, 28% had 3 shifts

#Here is the mean phylorate plot, which is for 3 shifts
plot.bammdata(ev1, lwd = 2)
addBAMMshifts(ev1, cex = 2) # The shift at the upper part of the tree corresponds to the Satyrinae and Charaxinae, the bottom shift is the Ithomiinae

# now the bayesian credible set of shifts
# read in the shift priors
pri1 <- read.table("prior_probs.txt", sep = ",", header = TRUE)

bayes_factor <- computeBayesFactors(mcmc1, pri1, burnin = 0.1) # 3 shifts has the highest Bayes factor (this is a table in the supplemental information) but there is no strong evidence for one model over another 

# plot a tree with the marginal shift probs, where branch length relates to the shift prob
Nym_marg_probs <- marginalShiftProbsTree(ev1)
plot.phylo(Nym_marg_probs, cex = 0.2) # this makes sense, because the lognest branchs correspond to the location of the rate shifts

# get credible sets of shift configurations
prior_shifts <- getBranchShiftPriors(Nym.pruned$phy, pri1)
cred_set <- credibleShiftSet(ev1, prior_shifts, set.limit = 0.95)
cred_set$number.distinct # 21 distinct rate shift configurations
summary(cred_set) # the majority of the credible sets all have 2 rate shifts
plot.credibleshiftset(cred_set) # all of the rate shifts plotted have shifts in the same general locations!

#now we find the "best" shift configuration
bestie <- getBestShiftConfiguration(ev1, prior = prior_shifts, BFcriterion = 5) 
par(mfrow = c(1, 1))
plot.bammdata(bestie, lwd = 2)
addBAMMshifts(bestie, cex = 2.5) # the upper node is the Satyrinae plus Charaxinae, the lower node is the Ithomiinae

# Bayes factors as evidence for rate shifts
Nym_bf <- bayesFactorBranches(ev1, prior_shifts)
plot.phylo(Nym_bf, cex = 0.2) # Now branch lengths are scaled to Bayes factors, with shorter branches having higher Bayes factors.
# nodelabels(format(round(Nym_bf$edge.length, 2), nsmall = 2), frame = "n", cex = 0.5, font = 2, adj = c(1.25, -1))

# as an alternative to the "best" tree we can find the tree with the maximum shift credibility
maxie <- maximumShiftCredibility(ev1, maximize = "product") # we will maximize the product of the branch-specific shifts
maxie_config <- subsetEventData(ev1, index = maxie$sampleindex)
plot.bammdata(maxie_config, lwd = 2)
addBAMMshifts(maxie_config, cex = 2) #this configuarion has rate shifts at the Satyrinae and the Itomiinae 

# Here is get an idea of the distribution of rates in the tree
Nym.rates <- getCladeRates(ev1)

par(mfrow=c(1, 2))
# plot speciation rate
hist(Nym.rates$lambda, las = 1, col = "grey", main = "", xlab=expression(paste(lambda)), ylim = c(0, 275), xlim = c(0.03, 0.055)) # looking at a histogram of speciation rates
abline(v=mean(Nym.rates$lambda), col = "red", lwd = 2) # mean
qlam <- quantile(Nym.rates$lambda, c(0.025, 0.975)) # 95% quantile CI
abline(v = qlam, col = "blue", lwd = 2)

#plot extinction rate
hist(Nym.rates$mu, las = 1, col = "grey", main = "", xlab = expression(paste(mu)), ylim = c(0, 275), xlim = c(0, 0.015)) # looking at a histogram of extinction rates
abline(v = mean(Nym.rates$mu), col = "red", lwd = 2) # mean
mean(Nym.rates$mu) # compute mean extinction rate
qmu <- quantile(Nym.rates$mu, c(0.025, 0.975)) # 95% quantile for mean
abline(v = qmu, col = "blue", lwd = 2)


#plot diversification rate
hist(Nym.rates$lambda - Nym.rates$mu, las = 1, col = "grey", main = "", xlab=expression(paste(epsilon)), ylim = c(0, 275), xlim = c(0.03, 0.053)) # looking at a histogram of speciation rates
abline(v=mean(Nym.rates$lambda - Nym.rates$mu), col = "red", lwd = 2) # mean
qlam <- quantile(Nym.rates$lambda - Nym.rates$mu, c(0.025, 0.975)) # 95% quantile CI
abline(v = qlam, col = "blue", lwd = 2)


# just for fun plot speciation and extinction estimates
plot(Nym.rates$lambda ~ Nym.rates$mu, pch = 19, las = 1, xlab = expression(paste(mu)), ylab = expression(paste(lambda)), ylim = c(0.035, 0.055), xlim = c(0, 0.015))
abline(lm(Nym.rates$lambda ~ Nym.rates$mu), lwd = 2) # note how low extinction rate estimates are relative to speciation rates. 
par(mfrow = c(1, 1))


Nym.dsc1 <- distinctShiftConfigurations(ev1, prior = Nym.prior.threshold, BFcriterion = 5)
plot.bammshifts(Nym.dsc1, ev1, rank = 1, legend = TRUE) # plot the best shift configuration onto the pretty tree, which has the satyrs and the itomiines

#Two last things to look at the is the speciation rate through time
Nym.mat <- getRateThroughTimeMatrix(ev1)

lam_Nym <- apply(Nym.mat$lambda, 2, quantile, seq(0.1, 0.9, by = 0.1))
vec1 <- max(Nym.mat$times) - Nym.mat$times
cols <- c('gray90', 'gray75', 'gray60', 'gray45')

plot.new()
par(mar = c(6 , 6, 1, 1))
plot.window(xlim = c(95, 0), ylim = c(0, 0.30))
for (i in 1:4){
	xcoords <- c(vec1, rev(vec1))
	ycoords <- c(lam_Nym[i, ], rev(lam_Nym[10-i, ]))
	polygon(x = xcoords, y = ycoords, col = cols[i], border = FALSE)
}
lines(x = vec1, y = colMeans(lam_Nym), lwd = 3, col = 'black')
axis(1, at = seq(95, -5, by = -5), cex.axis = 1.3)
axis(2, at = seq(0, 0.4, by=0.05), cex.axis=1.3, las=1)
mtext(side = 1, text = "Time before present (mya)", line = 3.5, cex = 1.2)
mtext(side = 2, text=expression(paste(lambda)), line=3.5, cex=1.2)


# and a pretty shift plot
library(gplots) # v2.14.1
pal <- rich.colors(10)
cexmin <- 2
cexmax <- 6
Nym.nodes <- getShiftNodesFromIndex(ev1, index = maxie$sampleindex)
Nym.cst <- cumulativeShiftProbsTree(ev1)
probvec <- numeric(length(Nym.nodes))

for (i in 1:length(Nym.nodes)){
	probvec[i] <- Nym_marg_probs$edge.length[Nym_marg_probs$edge[, 2] == Nym.nodes[i]]
}
size <- cexmin + probvec * (cexmax - cexmin)
colvec <- rep(pal[length(pal)], length(Nym.nodes))

cuts <- seq(0, 1, length.out = length(pal) + 1)

for (i in 1:(length(cuts) - 1)){
	for (k in 1:length(probvec)){
		if (probvec[k] > cuts[i] & probvec[k] <= cuts[i+1]){
			colvec[k] <- pal[i];
		}
	}	
}
lmat <- matrix(c(1, 1, 1, 1, 1, 1, 1, 2), nrow = 1)
edgewid <- rep(0.5, length(Nym.pruned$phy$edge.length))
edgecol <- rep('gray60', length(Nym.pruned$phy$edge.length))

edgewid[Nym.cst$edge.length > 0.9] <- 0.7
edgecol[Nym.cst$edge.length > 0.9] <- pal[10]
 
plot.new()
layout(lmat)
 
plot.phylo(Nym.pruned$phy, show.tip.label = FALSE, edge.width = 0.8, edge.color = edgecol, no.margin = TRUE, type = 'p', direction = 'upwards')
nodelabels(node = Nym.nodes, cex = size, pch = 21, bg = colvec)
 
plot.new()
par(mar=c(3,0,3,6))
plot.window(xlim = c(0, 4), ylim=c(-0.3, 1.3))
for (i in 1:(length(cuts)-1)){
	xco <- c(0, 2, 2, 0)
	yco <- c(cuts[i], cuts[i], cuts[i + 1], cuts[i + 1])
	polygon(x = xco, y = yco, col = pal[i])
}
axis(side=4,at=seq(0,1, by=0.2), cex.axis=2, font=2, las=1, pos=2.2)
mtext(text='Shift probability', side=4, line=4, cex=1) 


# looking at evolutionary cohort matrix, mostly because this is a pretty way of looking at the clades with increases in diversification rate
cohort_matrix <- getCohortMatrix(ev1)
cohorts(cohort_matrix, ev1, lwd = 3, pal = "temperature", use.plot.bammdata = TRUE, useraster = FALSE)


# getting mean clade-specific evolutionary rates to compare the clades on either side of the BAMM identified shift. Note this function uses the node as numbered by ape (find them using the nodelabels function)
#Functions to create 95% CIs
se <- function(x){
sqrt(var(x, na.rm = TRUE) / length(na.omit(x)))
} 

ci95 <- function(x){
	t.value <- qt(0.975, length(x) - 1)
	standard.error <- se(x)
	ci <- t.value * standard.error
	cat("mean" = mean(x), "95% CI =", mean(x) - ci, "to", mean(x) + ci, "\n")
}

# The "node 1 rate shift identified by BAMM corresponds to the Ithomiinae, which is node 386; it's sister node is 421
ithom_rates <- getCladeRates(ev1, node = 386, nodetype = "include") #Lots of the Ithomiini
mechan_rates <- getCladeRates(ev1, node = 421, nodetype = "include") #Mechanitina sensu (Brower et al. 2006)

ci95(ithom_rates$lambda)
quantile(ithom_rates$lambda, c(0.05, 0.95))

ci95(mechan_rates$lambda)
quantile(mechan_rates$lambda, c(0.05, 0.95))


t.test(ithom_rates$lambda, mechan_rates$lambda)

hist(ithom_rates$lambda, col = "#D3D3D37D", xlim = c(0.025, 0.16), ylim = c(0, 250), xlab = expression(paste(lambda)), main = "", las = 1, breaks = 20)
hist(mechan_rates$lambda, col = "#A9A9A97D", breaks = 20, add = TRUE)
legend("topright", legend = c(" Satyrini"," other Satyrinae"), pch = 15, col = c("#A9A9A97D", "#D3D3D37D"), bty = "n", pt.cex = 2)

ci95(ithom_rates$mu)
quantile(ithom_rates$mu, c(0.05, 0.95))

ci95(mechan_rates$mu)


# Part of Satyrini is node 595, sister node (other Satyrini) is node 716
satyr_rates <- getCladeRates(ev1, node = 595, nodetype = "include")
otherSatyrs <- getCladeRates(ev1, node = 716, nodetype = "include")

ci95(satyr_rates$lambda)
quantile(satyr_rates$lambda, c(0.05, 0.95)) # this is the 90% HPD

ci95(otherSatyrs$lambda) 
quantile(charax_rates$lambda, c(0.05, 0.95))

t.test(satyr_rates$lambda, otherSatyrs$lambda) # not sig different, but histograms overlap

# Function I found on by Nick Sabbe on StackOverflow to generate color codes based on names and desired transparency
makeTransparent<-function(someColor, alpha=100){
   newColor<-col2rgb(someColor)
   apply(newColor, 2, function(curcoldata){rgb(red = curcoldata[1],  green = curcoldata[2], blue=curcoldata[3], alpha = alpha, maxColorValue=255)})
}
makeTransparent("light grey", alpha = 125)
makeTransparent("dark grey", alpha = 125)

hist(otherSatyrs$lambda, col = "#D3D3D37D", xlim = c(0.025, 0.075), ylim = c(0, 250), xlab = expression(paste(lambda)), main = "", las = 1, breaks = 20)
hist(satyr_rates$lambda, col = "#A9A9A97D", breaks = 20, add = TRUE)
legend("topright", legend = c(" Satyrini"," other Satyrinae"), pch = 15, col = c("#A9A9A97D", "#D3D3D37D"), bty = "n", pt.cex = 2)

ci95(satyr_rates$mu)
quantile(satyr_rates$mu, c(0.05, 0.95)) 

ci95(otherSatyrs$mu)
quantile(otherSatyrs$mu, c(0.05, 0.95))

t.test(satyr_rates$mu, otherSatyrs$mu) # Big difference

hist(otherSatyrs$mu, col = "#D3D3D37D", xlim = c(0, 0.032), ylim = c(0, 400), xlab = expression(paste(mu)), main = "", las = 1, breaks = 20)
hist(satyr_rates$mu, col = "#A9A9A97D", breaks = 40, add = TRUE)
legend("topright", legend = c(" Satyrini"," other Satyrinae"), pch = 15, col = c("#A9A9A97D", "#D3D3D37D"), bty = "n", pt.cex = 2)

mean(satyr_rates$lambda - satyr_rates$mu) # Diversification rate
mean(otherSatyrs$lambda - otherSatyrs$mu) # is the explanation that the remaining Satyrinae have a reduced extinction rate here? And that it is statistically significant (but the biological significance is unknown)






# Here we will comapre the host breadth of clades with increased diversification rates and their sister clades with lower rates. First we need to extract the subtrees

idNym <- id.subtrees(Nym.pruned$phy)
# names(idNym)
# pdf(file='Nym_nodes_Hamm.pdf', bg = "white")
# plot(idNym$subtree[[1]],show.node.label = TRUE, cex = 0.2) This plot has all nodes numbered
# dev.off()
Nymsubtree <- idNym$subtree


# PIC and Wilcoxan sign rank test identified from turboMEDUSA run tMED1 and BAMM!

# split 2 is at node 26, it's sister is node 20, and the parent node is 19. Node 26 has a speed up relative to the background rate.
Nym.19 <- treedata(phy = Nymsubtree[[19]], data = Covs, sort = TRUE, warnings = FALSE)
pic.19 <- pic(x = Nym.19$data[, 1], phy = Nym.19$phy, scaled = TRUE, var.contrasts = FALSE)
#plot(Nym.19$phy) ,nodelabels(round(pic.19, digits = 4), frame = "n")
# now we ask of the node with the higher diversification rate has more hosts than it sister, pic.19[2:7] is the clade with slower div rate. Note that in all comparisons we extract the contrast value for the branch that connect the two clades
hist(pic.19[8:21], col = "grey", xlim = c(-0.5, 0.5), main = "Node 2", , xlab = "Corrected host breadth", las = 1)
hist(pic.19[2:7], col = "green", add = TRUE)
legend("topright", legend  = c("Fast clade", "Slow clade"), pch = 15, col = c("grey", "green"), bty = "n", cex = 1.5)

library(exactRankTests) # v0.8-27
wilcox.exact(pic.19[8:21], pic.19[2:7], exact = TRUE, alternative = "less") # specifically asking if the "faster" clade has higher host breadth than the "slower" clade
wilcox.exact(pic.19[8:21], pic.19[2:7], exact = TRUE, alternative  =  "two.sided") # This gives no evidence against identical distributions


# split 3 is at node 217 (pic.216[2:121]), its sister node is 338, parent node is 216
Nym.216 <- treedata(phy = Nymsubtree[[216]], data = Covs, sort = TRUE, warnings = FALSE)
pic.216 <- pic(x = Nym.216$data[,1], phy = Nym.216$phy, scaled = TRUE, var.contrasts = FALSE)
# plot(Nym.216$phy); nodelabels(round(pic.216, digits = 4), frame = "n")
wilcox.exact(pic.216[2:122], pic.216[123:162], exact = TRUE, alternative = "less")
wilcox.exact(pic.216[2:121], pic.216[122:161], exact = TRUE, alternative = "two.sided")

# split 4 is on node 138, its sister is node 129, parent node is 128
# pic.128[2:10] is the faster clade
Nym.128 <- treedata(phy=Nymsubtree[[128]], data = Covs, sort = TRUE, warnings = FALSE)
pic.128 <- pic(x=Nym.128$data[, 1], phy = Nym.128$phy, scaled = TRUE, var.contrasts = FALSE)
wilcox.exact(pic.128[2:10], pic.128[10:12], exact = TRUE, alternative = "less")
# plot(Nym.128$phy); nodelabels(round(pic.128, digits = 4), frame="n")
wilcox.exact(pic.128[2:10], pic.128[10:12], exact = TRUE, alternative = "two.sided")


# split 5 is on node 236 (a slowdown from background rate), its sister is 231, parent node is 230
# pic.230[6:8] is the fast clade
Nym.230 <- treedata(phy = Nymsubtree[[230]], data = Covs, sort = TRUE, warnings = FALSE)
pic.230 <- pic(x=Nym.230$data[, 1], phy=Nym.230$phy, scaled = TRUE, var.contrasts = FALSE)
# plot(Nym.230$phy); nodelabels(round(pic.230, digits = 4), frame = "n")
wilcox.exact(pic.230[6:8], pic.230[2:5], exact = TRUE, alternative = "less")
wilcox.exact(pic.230[6:8], pic.230[2:5], exact = TRUE, alternative = "two.sided")

# Previous studies have reported a correlation between butterfly taxon size and host breadth. Let's make sure we can recreat that finding. PIC and then look at correlation between host breadth and taxon size.
pic.Host <- pic(x = Hosts, phy = Nym.pruned$phy, scaled = TRUE, var.contrasts=FALSE)
pic.Gen <- pic(x = Nsp, phy = Nym.pruned$phy, scaled = TRUE, var.contrasts = FALSE)

super.pic <- lm(pic.Host ~ pic.Gen) # pic is the way to go, makes the data nomrally distributed
summary(super.pic) # small slope and small R-squared (0.14). Statistically significant, however
par(mfrow = c(2, 2))
plot(super.pic) # residuals look decent, no shotgun pattern of resids v fitted values
par(mfrow = c(1, 1))

# Legendre and Desdevises recommend a regression through the origin and a premutation test to correct for the bias

(pic.origin <- lmorigin(pic.Host ~ pic.Gen, nperm = 10000)) # an R^2 of 0.14, no different than a normal lm


# always a good idea to visualize your data
hist(pic.Host, col = "grey", las = 1)
hist(pic.Gen, col = "grey", las = 1)

cor.test(pic.Host, pic.Gen, method = "spearman") # For supporting material, Spearman's rho is 0.4

library(arm) #v 1.6-10 #just a simple function to plot a 95%CI based on resampling a regression. I got this idea from Ian Dworkin at Michigan State. If he tells me he read this I'll buy him a beer for the idea. 
simmie <- function(){
	sim.lm <- sim(super.pic, n.sims = 1)
	abline(a = coef(sim.lm)[, 1], b = coef(sim.lm)[, 2], col = "dark grey", lty = 3, lwd = 1)
}

# pdf(file = "pic_1a.pdf", bg = "white")
plot(pic.Gen, pic.Host, pch = 19, las = 1, ylab = "Host breadth", xlab = "Genus size")
invisible(replicate(95, simmie())) # 95% CIs
abline(super.pic, lwd = 3, lty = 1)
# dev.off()





# here we look at specialists and generalists in terms of clades, rather than tips. Read in file with butterfly genera and the host plant family that each genus feeds on (one row per host record)
library(reshape2) # v1.04
Nym_host <- read.csv("Nym_host_9-Sept-2014.csv", header = TRUE)
Nym_host1 <- melt(Nym_host, id.var = "Bfly")
dat <- as.data.frame.matrix(with(Nym_host1, table(Bfly, value)))
dim(dat)

# There are two additional genera in these data compare the the "Covs" file we loaded at the beinging of this analysis, otherwise they containt the same genera of butterflies. 
setdiff(row.names(dat), row.names(Covs))

library(vegetarian) # v1.2
datt <- dat
datt$X <- row.names(dat)
subs <- id.subtrees(Nym.pruned$phy)

# the rate shifts have been identified at the following locations:
# Node 2: corresponds to subtree #26 (faster rate), sister subtree is #20 
# Node 3: subtree 217 (faster*), sister subtree is 231
# Node 4: subtree 138 , sister subtree is 129 (faster)
# Node 5: subtree 236 , sister is 231 (faster)

# This is a function that we wrote to calculate Jost's D for each clade and derive an uncorrected P value asking if they two host breadths are equivalent.

phylo.beta.d <- function(fast.tree, slow.tree, N){
	intree <- match(subs$subtree[[fast.tree]]$tip.label, datt$X)
	outree <- match(subs$subtree[[slow.tree]]$tip.label, datt$X)
	dat.in <- dat[intree, ]
	dat.out <- dat[outree, ]
	dat.in <- dat.in[, -which(apply(dat.in, 2,sum) == 0)]
	dat.out <- dat.out[, -which(apply(dat.out, 2, sum) == 0)]
	print(d(dat.in, lev = "beta", q = 0, boot = TRUE, boot.arg = list(num.iter = N)))
    print(d(dat.out, lev = "beta", q = 0, boot = TRUE, boot.arg = list(num.iter = N)))
	samps <- min(dim(dat.in)[1], dim(dat.out)[1])
	in.boot <- NA; out.boot <- NA
		for(i in 1:N){
			in.boot[i] <- d(dat.in[sample(1:dim(dat.in)[1], samps, replace=TRUE), ], 	lev = "beta", q = 0)
			out.boot[i] <- d(dat.out[sample(1:dim(dat.out)[1], samps, replace = 		TRUE), ], lev = "beta", q = 0)	
	}
		the.faster <- in.boot
		the.slower <- out.boot
		(sum(the.slower - the.faster < 0) / N)
}

phylo.beta.d(26, 20, 10000)
phylo.beta.d(217, 338, 10000)
phylo.beta.d(129, 138, 10000)
phylo.beta.d(231, 236, 10000)

p.adjust(c(0.5215, 0.0402, 0.8427, 0.6701), method = "holm") #order of p-values entered = node 1, node 2, node 3, node 4, node 5

# As a fun exercise, it is cool to see how host use actually seems fairly conserved within a clade. Let's look at the Ithomiines from Node 2 (branch 26). 

Ithom <-  match(subs$subtree[[217]]$tip.label, datt$X)
data.Ithom <- dat[Ithom, ] #This generates a matrix with EVERY host plany
data.Ithom <- data.Ithom[, -which(apply(data.Ithom, 2, sum) == 0)] #This removes host plants not used by the clade in question.
data.Ithom # Cool! See how nearly everything in that clade feeds on Solanaceae! Try this with the satyrs (branch 217) and confirm that nearly every genus in that clade feeds on Poaceae






###################################################
# now looking at phylogenetic signal as it relates to host breadth

# Bin the data into binary sets based on the number of host plant families the butterfly genus feeds on
# specialists feed on ONE plant family, generalists on > 1
# convert to 0's and 1's
W <- Hosts
head(W)

bi=NA
for(t in 1:length(W)){
	if (W[t] == 1) bi[t] = 0
	else bi[t] = 1
}
names(bi) <- row.names(Covs)
head(bi)
sum(bi) # 173 polyphagous genera (206 specialist)

# specialists feed on TWO or fewer host families, generalists on > 2
tri = NA
for(t in 1:length(W)){
	if (W[t] <= 2) tri[t] = 0
	else tri[t] = 1
}
names(tri) <- row.names(Covs)
head(tri)
sum(tri) # 106 generalist genera (271 specialist)

# specialist feed on THREE or fewer host plant families, generalists on > 3
quad = NA
for(t in 1:length(W)){
	if (W[t] <= 3) quad[t] = 0
	else quad[t] = 1
}
names(quad) <- row.names(Covs)
head(quad)
sum(quad) # 70 generalist genera, 307 specialist genera

# specialists feed on FOUR or fewer host plant families, generalists on > 4
quint = NA
for(t in 1:length(W)){
	if (W[t] <= 6) quint[t] = 0
	else quint[t] = 1
}
names(quint) <- row.names(Covs)
head(quint)
sum(quint) # 30 generalist genera, 347 specialist taxa


# calculating phylogenetic signal based on data partitions, here using Blomberg's K with simulations to estimate the parameters
(Host.K <- phylosig(Nym.pruned$phy, Hosts, method = "K", test = TRUE, nsim = 10000)) #This method is for continuous characters but it is fine for host plants
(K.bi <- phylosig(Nym.pruned$phy, bi, method = "K", test = TRUE, nsim = 10000))
(K.tri <- phylosig(Nym.pruned$phy, tri, method = "K", test = TRUE, nsim = 10000))
(K.quad <- phylosig(Nym.pruned$phy, quad, method = "K", test = TRUE, nsim = 10000))
(K.quint <- phylosig(Nym.pruned$phy, quint, method = "K", test = TRUE, nsim = 10000))


# Pagel's lambda
(Host.lambda <- phylosig(Nym.pruned$phy, Hosts, method = "lambda", test = TRUE))
(lambda.bi <- phylosig(Nym.pruned$phy, bi, method = "lambda", test = TRUE))
(lambda.tri <- phylosig(Nym.pruned$phy, tri, method = "lambda", test = TRUE))
(lambda.quad <- phylosig(Nym.pruned$phy, quad, method = "lambda", test = TRUE))
(lambda.quint <- phylosig(Nym.pruned$phy, quint, method = "lambda", test = TRUE))

# more phylogenetic signal fun. This time with the OU model and our power to reject a pure BM model based on the above data partitions
OU_cont <- fitContinuous(Nym.pruned$phy, dat = Nym.pruned$data[, 1], model = "OU")
OU_bi <- fitContinuous(Nym.pruned$phy, dat = bi, model = "OU")
OU_tri <- fitContinuous(Nym.pruned$phy, dat = tri, model = "OU")
OU_quad <- fitContinuous(Nym.pruned$phy, dat = quad, model = "OU")
OU_quint <- fitContinuous(Nym.pruned$phy, dat = quint, model = "OU")

# now run the power analysis 
# library(pmc) # v0.0-8
# library(snowfall) #v 1.84-6
# sfInit(parallel = TRUE, cpus = 6)
# tree_pow_1 <- treepower(ape2ouch(Nym.pruned$phy), nboot = 1000, threshold = 0.95, alpha = c(0.02, 0.02459, 0.03 ), cpu = 6); sfStop()

# sfInit(parallel = TRUE, cpus = 6)
# tree_pow_2 <- treepower(ape2ouch(Nym.pruned$phy), nboot = 1000, threshold = 0.95, alpha = c(0.02, 0.03036498, 0.04 ), cpu = 6); sfStop()

# sfInit(parallel = TRUE, cpus = 6)
# tree_pow_3 <- treepower(ape2ouch(Nym.pruned$phy), nboot = 1000, threshold = 0.95, alpha = c(0.02, 0.02488918, 0.03 ), cpu = 6); sfStop()

# sfInit(parallel = TRUE, cpus = 6)
# tree_pow_4 <- treepower(ape2ouch(Nym.pruned$phy), nboot = 1000, threshold = 0.95, alpha = c(0.02, 0.02655768, 0.03 ), cpu = 6); sfStop()

# sfInit(parallel = TRUE, cpus = 6)
# tree_pow_5 <- treepower(ape2ouch(Nym.pruned$phy), nboot = 1000, threshold = 0.95, alpha = c(0.02, 0.03659707, 0.04 ), cpu = 6); sfStop()

###################################################
# examining the effects of host breadth on rates
library(diversitree) # v0.9-6
# BiSSE analysis, based on data bins, in both MLE and Bayesian frameworks

# for "bi" data set
# we have 378/500 (0.756), is (0.244 missing taxa) genera represented in our tree. The proportion of missing taxa for each character states should correlate with the proportion of taxa with that character state. 

# for the "bi" state, 173 generalist
length(bi) - sum(bi) #205 specialists

173/378 # proportion of generalist taxa
205/378 # proportion of specialist taxa
# this assumes sampling fraction is dependent on character state

1 - (0.244) * (173 / 378) # I think we have ~88% of the specialists
1 - (0.244) * (205 / 378) # I think we have ~86% of the generalists

samp.bi <- as.numeric(c(1 - (0.244) * (173 / 378), 1 - (0.244) * (205 / 378)))
samp.bi


# using sampling.f to accont for missing taxa
bisse.bi.1a <- make.bisse(tree = Nym.pruned$phy, states = bi, sampling.f = samp.bi)

start.bi.1a <- starting.point.bisse(tree = Nym.pruned$phy, yule = FALSE, q.div = 0.6) #yule=FALSE because we do not want a pure birth model
fit.bi.1a <- find.mle(bisse.bi.1a, start.bi.1a, method = "subplex")

prior5 <- make.prior.uniform(lower = 0, upper = 10) # uniform prior!
# tmp_2a <- mcmc(bisse.bi.1a, fit.bi.1a$par, nsteps = 100, w = 0.1)
w_2a <- diff(sapply(tmp_2a[2:7], quantile, c(0.05, 0.95))) # get values for w

# bisse.mcmc.bi.1a <- mcmc(bisse.bi.1a, fit.bi.1a$par, nsteps = 10000, prior = prior5, lower = 0, w = w_2a, print.every = 100)#This took about three hours


# for "tri" data set
length(tri) - sum(tri) #272 specialists, 106 generalists
1 - (0.244) * (271 / 378) # I think we have 82% of the specialists
1 - (0.244) * (107 / 378) # I think we have 93% of the specialists
samp.tri <- as.numeric(c(1 - (0.244) * (272 / 378), 1 - (0.244) * (106 / 378))) # note we am trying to estimate the proportion of missing taxa uniquely for each bin

bisse.tri.1a <- make.bisse(tree = Nym.pruned$phy, states = tri, sampling.f = samp.tri)
start.tri.1a <- starting.point.bisse(tree = Nym.pruned$phy, q.div = 0.6)
fit.tri.1a <- find.mle(bisse.tri.1a, start.tri.1a, method = "subplex")

# tmp_3a <- mcmc(bisse.tri.1a, fit.tri.1a$par, nsteps = 100, w = 0.1)
w_3a <- diff(sapply(tmp_3a[2:7], quantile, c(0.05, 0.95)))

# bisse.mcmc.tri.1a <- mcmc(bisse.tri.1a, fit.tri.1a$par, nsteps = 10000, prior = prior5, lower = 0, w = w_3a, print.every = 100)


# for "quad" data set
length(quad) - sum(quad) # 70 generalist taxa, 308 specialists
1 - (0.244) * (308 / 378) # 80% of the specialists
1 - (0.244) * (70 / 378) # 95% of the generalists

samp.quad <- as.numeric(c(1 - (0.244) * (308 / 378), 1 - (0.244) * (70 / 378)))

bisse.quad.1a <- make.bisse(tree = Nym.pruned$phy, states = quad, sampling.f = samp.quad)
start.quad.1a <- starting.point.bisse(tree = Nym.pruned$phy, q.div = 0.5) # No matter what I did here with q.div I was not able to ge the MLE close to what the mcmc came up with, so I ended up feeding it unique starting values to come up with an MLE that was within the MCMC posterior distribution. 
start.quad.1b <- c(0.001, 0.008, 0.0, 0.0, 0.11577877, 0.11577877)
names(start.quad.1b) <- names(start.quad.1a)

fit.quad.1a <- find.mle(bisse.quad.1a, start.quad.1b, method = "subplex") # something a bit off here, misses the mle. Really work this. 

# tmp_4a <- mcmc(bisse.quad.1a, fit.quad.1a$par, nsteps = 100, w = 0.1)
w_4a <- diff(sapply(tmp_4a[2:7], quantile, c(0.05, 0.95)))

# bisse.mcmc.quad.1a <- mcmc(bisse.quad.1a, fit.quad.1a$par, nsteps = 10000, prior = prior5, lower = 0, w = w_4a, print.every = 100)


# for "quint" data set
length(quint) - sum(quint) # 348 specialists, 30 generalists (must have lots of taxa)
1 - (0.244) * (348 / 378) # 77% of the specialists
1 - (0.244) * (30 / 378) # 98% of the generalists

samp.quint <- as.numeric(c(1 - (0.244) * (348 / 378), 1 - (0.244) * (30 / 378)))

bisse.quint.1 <- make.bisse(tree = Nym.pruned$phy, states = quint, sampling.f = samp.quint)
start.quint.1 <- starting.point.bisse(tree = Nym.pruned$phy, q.div = 0.3)
fit.quint.1 <- find.mle(bisse.quint.1, start.quint.1, method = "subplex")

# tmp_5 <- mcmc(bisse.quint.1, fit.quint.1$par, nsteps = 100, w = 0.1)
w_5 <- diff(sapply(tmp_5[2:7], quantile, c(0.05, 0.95)))

# bisse.mcmc.quint.1 <- mcmc(bisse.quint.1, fit.quint.1$par, nsteps = 10000, prior = prior5, lower = 0, w = w_5, print.every = 100)

# it looks like the model with variable rates for lambda and inertia is best.
# let's look at correlations
library(adephylo) # v1.1-6
library(phylobase) # v0.6.5.2
Nym.phylo4d <- phylo4d(Nym.pruned$phy, data.frame(Hosts, Nsp))
abouheif.moran(Nym.phylo4d) # looks like a strong correlation 

# Let's look at a test of phylogenetic inertia based on Ollier et al. (2005)
ornor <- as.matrix(orthobasis.phylo(Nym.pruned$phy))#compute the orthonormal basis for a phylogeny using phylogenetic structure
two.ornor <- ornor[, 1:2] # must use first two columns

ornor.lm <- anova(lm(Nym.pruned$data[, 1] ~ two.ornor)) 
ornor.lm # look for inertia in hosts. Yup, it looks like it is there
ortho.1 <- adephylo::orthogram(Nym.pruned$data[, 1], Nym.pruned$phy, nrepet = 10000)
# nodes with var decomp > 0.02: 3, 5, 12, 75, 83, 106, 122, 134, 206, 209, 288, 298 - to get these node I cut out the orthogram code and ran it line by line, then called which(w$phylogram > 0.02).

# for which(w$phylogram > 0.03)
# 3   5  12  83 206 209 298 

# What is the diet breadth (in terms of host plant families) for these nodes? We just want Jost's d
# For node 5
intree5 <- match(subs$subtree[[5]]$tip.label, datt$X); length(intree5)
dat.in5 <- dat[intree5, ]
dim(dat.in5)
d5 <- d(dat.in5, lev = "beta", q = 1, boot = TRUE, boot.arg = list(num.iter = 1000)) # d = 2.13, se = 0.18

# function to permute the data and come up with a null distribution
permuD <- function(dat.in, iterations){
	tdat.in <- dat.in
	tdat.in <- tdat.in[, -c(which(apply(tdat.in, 2, sum) == 0))]
	nullD <- NA
	for(iter in 1:iterations){
 		if (iter %% 100 == 0){cat("\n", iter, "of",iterations)}
	temp <- tdat.in
	for(i in 1:dim(tdat.in)[1]){
		temp[i, ] <- sample(tdat.in[i, ], size = dim(tdat.in)[2])
		}
	nullD[iter] <- d(temp, lev = "beta", q = 1)
	}
	return(nullD)
}

test5 <- permuD(dat.in5, iterations = 1000)

# what is the mean host breadth for this group?
Nym.sub5 <- treedata(phy = Nymsubtree[[5]], data = Covs, sort = TRUE, warnings = FALSE)
mean(Nym.sub5$data[, 1]) # mean 1.26,
# std <- function(x) sd(x)/sqrt(length(x)) #simple function to calculate the se

# For node 83
# Nym.sub83 <- treedata(phy=Nymsubtree[[83]], data = Covs, sort = TRUE, warnings = TRUE)
intree83 <- match(subs$subtree[[83]]$tip.label, datt$X); length(intree83)
dat.in83 <- dat[intree83, ]
d83 <- d(dat.in83, lev = "beta", q = 1, boot = TRUE, boot.arg = list(num.iter = 1000)) # d=9.7, se = 0.95

test83 <- permuD(dat.in83, iterations = 1000)


# For node 206
intree206 <- match(subs$subtree[[206]]$tip.label, datt$X); length(intree206)
dat.in206 <- dat[intree206, ]
d206 <- d(dat.in206, lev = "beta",q = 1, boot = TRUE, boot.arg = list(num.iter = 1000)) # d = 2.34, se = 0.43

test206 <- permuD(dat.in206, iterations = 1000)

Nym.sub206 <- treedata(phy=Nymsubtree[[206]], data = Covs, sort = TRUE, warnings = FALSE)

# For node 253
intree253 <- match(subs$subtree[[253]]$tip.label, datt$X); length(intree253)
dat.in253 <- dat[intree253, ]
d253 <- d(dat.in253,lev = "beta", q = 1, boot = TRUE, boot.arg = list(num.iter = 1000)) # d = 1.51 se = 0.10

test253 <- permuD(dat.in253, iterations = 1000)

Nym.sub253 <- treedata(phy = Nymsubtree[[253]], data = Covs, sort = TRUE, warnings = FALSE)

# for node 122
intree122 <- match(subs$subtree[[122]]$tip.label, datt$X); length(intree122)
dat.in122 <- dat[intree122, ]
d122 <- d(dat.in122, lev = "beta", q = 1, boot = TRUE, boot.arg = list(num.iter = 1000)) # mean = 2.89, se = 0.29

test122 <- permuD(dat.in122, iterations = 1000)

Nym.sub122 <- treedata(phy=Nymsubtree[[122]], data = Covs, sort = TRUE, warnings = FALSE)

# pdf(file="permu_test.pdf")
par(mfrow = c(5,1))
hist(test5, xlim = c(2, 10), col = "grey", las = 1, main = "Node 5", xlab = "d")
abline(v = d5$D.Value, col = "red", lwd = 2)

hist(test83, xlim = c(8, 14), col = "grey", las = 1, main = "Node 83", xlab = "d")
abline(v = d83$D.Value, col = "red", lwd = 2)

hist(test206, xlim = c(0, 5), col = "grey", las = 1, main = "Node 206", xlab = "d")
abline(v = d206$D.Value, col = "red", lwd = 2)

hist(test253, xlim = c(0, 14), col = "grey", las = 1, main = "Node 253", xlab = "d")
abline(v = d253$D.Value, col = "red", lwd = 2)

hist(test122, xlim = c(2, 4), col = "grey", las = 1, main = "Node 122", xlab = "d")
abline(v = d122$D.Value, col = "red", lwd = 2)
# dev.off()


# Now running the MuSSE, an extension of BiSSE that allows for multiple states. I'm setting the maximum for host breadth to 6 because k=23 will take forever. I can't think of anyone that would claim feeding on 6 host plants is not a generalist state.
state6 = NA
for(t in 1:length(Hosts)){
	if (Hosts[t] >= 6) state6[t] = 6
	else state6[t] = Hosts[t]
	}
names(state6) <- names(Hosts)



# looks like big transition from 1-6 
musse.1 <- make.musse(tree = Nym.pruned$phy, k = 6,states=state6, sampling.f = 0.756, strict = FALSE)
start.musse.1 <- starting.point.musse(Nym.pruned$phy, k = 6, yule = FALSE)

 fit.musse.1b <- find.mle(musse.1, start.musse.1, method = "subplex", control = list(maxit = 100000))




# making the figures in the paper
# pdf(file="Pub_fig_1b.pdf", bg="white")
par(mar = c(1,1,1,1), fg = "transparent")
plotTree(Nym.pruned$phy, lwd = 0.4, use.edge.length=TRUE, edge.width = 1.5)
segments(rep(97, ngenera), 1 : ngenera, rep(97, ngenera) + (0.4 * Hosts), 1:ngenera, lwd = 0.4, col = "black") # plot the tree
par(fg = "black")
text(97, -4, "Host plants", font = 2, cex = 0.5)
segments(rep(108, ngenera), 1:ngenera, rep(108, ngenera) + (0.04 * Nsp), 1:ngenera, lwd = 0.4, col = "black")
text(108, -4, "Genus size", font = 2, cex = 0.5)


# plot the iteRates identified switch points
# node 1 - this is the background rate for the tree
# points(0.1, 50, pch = 19, col = "black", cex = 1)
# text(-1.5, 55, "1)", cex = 0.7)
# node 2
points(72, 28, pch = 19, col = "black", cex = 1)
text(69, 28, "1)", cex = 0.7)
# node 3
points(38, 254.5, pch = 19, col = "black", cex = 1)
text(36, 260, "2)", pch = 19, col = "black", cex = 0.7)
# node 4
points(52.5, 135, pch = 19, col = "black", cex = 1)
text(49.5, 136, "3)", cex = 0.7)
# node 5
points(58, 233, pch = 19, col = "black", cex = 1)
text(54, 232, "4)", cex = 0.7)

# add the segments for the clades identified by orthonormal decomp.
segments(95.7, 340, 95.7, 250, lwd = 2, lend = "square")
segments(95.7, 209, 95.7, 199, lwd = 2, lend = "square")
segments(95.7, 116, 95.7, 123, lwd = 2, lend = "square")
segments(95.7, 80, 95.7, 104, lwd = 2, lend = "square")
segments(95.7, 2, 95.7, 44, lwd = 2, lend = "square")


# add the histograms
par(fig = c(0.05, 0.14, 0.84, 1), new = TRUE)
plot.new()
hist(Covs$Nhost, ylim = c(0, 300), xlim = c(0, 25), col = "grey", las = 1, xlab = "", main = '', breaks = 20, ylab = "", cex.axis = 0.5, mgp = c(0, 0.5, 0))
text(8, 250, labels = "A)")

par(fig = c(0.18, 0.27, 0.84, 1), new = TRUE)
hist(Covs$Nsp, ylim = c(0, 350),  xlim = c(0, 200), col = "grey", las = 1, xlab = "", main = "", ylab="", cex.axis = 0.5, mgp = c(0, 0.5, 0))
text(75, 300, labels = "B)")
# dev.off()


# Figure 2 (posterior plots from bisse mcmc) we found that looking at density plots was too confusing in grey scale so we went with the violin plots. We wrote these functions to create the plots. If you want to see these as density plots email CAH and I'll send you the code. 
library(vioplot) # v0.2

jafviol<-function(x, wex = 0.5, ltyCI = 2, MLE = NULL, ltyMLE = 1){
	# MLE should be in order bil0, bil1
	par(bty = "l")
	bir0 <- x$lambda0 - x$mu0
	bir1 <- x$lambda1 - x$mu1
	bil0 <- x$lambda0
	bil1 <- x$lambda1
	maxy <- max(c(bir0, bir1, bil0, bil1)) + 0.01
	miny <- min(c(bir0, bir1, bil0, bil1)) - 0.01
	plot(seq(from = miny, to = maxy, length = 7), seq(from = 0.5, to = 5.5, length = 7), type = "n", las = 1, yaxt = "n", ylab = "", xlab = "", xaxt = "n")
	axis(side = 1, padj = -0.5)

	vioplot(bir0, at = 1, wex = wex, add = TRUE, drawRect = FALSE, col = "white", horizontal = TRUE); segments(quantile(bir0, 0.975), 0.55, quantile(bir0, 0.975), 1.45 , lty = ltyCI); segments(quantile(bir0, 0.025), 0.55, quantile(bir0, 0.025), 1.45, lty = ltyCI)

	vioplot(bir1, at = 2, wex = wex, add = TRUE, drawRect = FALSE, col = "dark grey", horizontal = TRUE); segments(quantile(bir1, 0.975), 1.55, quantile(bir1, 0.975), 2.45, lty = ltyCI); segments(quantile(bir1, 0.025), 1.55, quantile(bir1, 0.025),2.45, lty = ltyCI)

	vioplot(bil0, at = 4, wex = wex, add = TRUE, drawRect = FALSE, col = "white", horizontal = TRUE)
	segments(quantile(bil0, 0.975), 3.55, quantile(bil0, 0.975), 4.45, lty = ltyCI)
	segments(quantile(bil0, 0.025), 3.55, quantile(bil0, 0.025), 4.45, lty = ltyCI)

	vioplot(bil1, at = 5, wex = wex, add = TRUE, drawRect = FALSE, col = "dark grey", horizontal = TRUE); segments(quantile(bil1, 0.975), 4.55, quantile(bil1, 0.975), 5.45, lty = ltyCI); segments(quantile(bil1, 0.025), 4.55, quantile(bil1, 0.025),5.45, lty = ltyCI)
	segments(MLE[1], 3.50,MLE[1], 4.5, lty = ltyMLE)
	segments(MLE[2], 4.50,MLE[2], 5.5, lty = ltyMLE)

	mtext(c(expression(paste("r"[0])), expression(paste("r"[1]))), at = c(1,2), side = 2, line = 0.5, cex = 1, las = 1)
	mtext(expression(paste(lambda[0])), side = 2, line = 0.7, cex = 1, at = 4, las = 1)
	mtext(expression(paste(lambda[1])), side = 2, line = 0.7, cex = 1, at = 5, las = 1)
	# mtext("parameter value",side = 1,cex = 0.7,line = 3)
}

jafviolq<-function(x, wex = 0.5, ltyCI = 2, MLE = NULL, ltyMLE = 1){
	par(bty = "l")
	q01 <- x$q01
	q10 <- x$q10
	maxy <- max(c(q01, q10)) + 0.01
	miny <- min(c(q01 ,q10)) - 0.01
	plot(seq(from = miny, to = maxy, length = 5), 0:4, type = "n", las = 1, yaxt = "n", ylab = "",xlab = "", xaxt = "n", bty = "n")
	axis(side = 1, padj = -0.5)
	vioplot(q01, at = 1, wex = wex, add = TRUE, drawRect = FALSE, col = "white", horizontal = TRUE) 
	segments(quantile(q01, 0.975), 0.55, quantile(q01, 0.975), 1.45, lty = ltyCI)
	segments(quantile(q01, 0.025), 0.55,quantile(q01, 0.025), 1.45, lty = ltyCI)
	segments(MLE[1], 0.50, MLE[1], 1.5, lty = ltyMLE)

	vioplot(q10, at = 3, wex = wex, add = TRUE, drawRect = FALSE, col = "dark grey", horizontal = TRUE); segments(quantile(q10, 0.975), 2.55, quantile(q10, 0.975),3.45, lty = ltyCI); segments(quantile(q10, 0.025), 2.55, quantile(q10, 0.025),3.45, lty = ltyCI)
	segments(MLE[2], 2.50, MLE[2], 3.5, lty = ltyMLE)
	mtext(c(expression(paste("q"["01"])), expression(paste("q"["10"]))), at = c(1,3), side = 2,line = 0.7, cex = 1, las = 1)
	# mtext("parameter value",side=1,cex=1,line=3)
}
	
# pdf(file = "violplot_JAF2.pdf", bg = "white")	
# quartz(width = 7,height = 6)
ltyCIs <- 3	
par(mar=c(2, 2, 0, 0))
par(fig=c(0.05, 0.45, 0.80, 1))
jafviol(bisse.mcmc.bi.1a, wex = 1.2, ltyCI = ltyCIs, MLE = mlez[1, c(1, 2)])
par(fig = c(0.55, 0.95, 0.80, 1), new = TRUE)
jafviolq(bisse.mcmc.bi.1a, wex = 1.2, ltyCI = ltyCIs, MLE = mlez[1, c(5, 6)])
par(fig = c(0.05, 0.45, 0.55, 0.75), new = TRUE)
jafviol(bisse.mcmc.tri.1a, wex = 1.2, ltyCI = ltyCIs, MLE = mlez[2, c(1, 2)])
par(fig = c(0.55, 0.95, 0.55, 0.75), new = TRUE)
jafviolq(bisse.mcmc.tri.1a, wex = 1.2, ltyCI = ltyCIs, MLE = mlez[2, c(5, 6)])
par(fig = c(0.05, 0.45, 0.30, 0.5), new = TRUE)
jafviol(bisse.mcmc.quad.1a,wex = 1.2, ltyCI = ltyCIs, MLE = mlez[3, c(1, 2)])
par(fig = c(0.55, 0.95, 0.30, 0.5), new = TRUE)
jafviolq(bisse.mcmc.quad.1a,wex = 1.2, ltyCI = ltyCIs, MLE = mlez[3, c(5, 6)])
par(fig = c(0.05, 0.45, 0.05, 0.25), new = TRUE)
jafviol(bisse.mcmc.quint.1, wex = 1.2, ltyCI = ltyCIs, MLE = mlez[4, c(1, 2)])
par(fig = c(0.55, 0.95, 0.05, 0.25), new = TRUE)
jafviolq(bisse.mcmc.quint.1, wex = 1.2, ltyCI = ltyCIs, MLE = mlez[4, c(5, 6)])
par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), new = TRUE)
plot(0:1, 0:1, xaxt = "n", yaxt = "n", ylab = "", xlab = "", type = "n")
text(0.25, 0, "Parameter value", cex = 1.5)
text(0.8, 0, "Parameter value", cex = 1.5)
ledders <- c("A", "B", "C", "D", "E", "F", "G", "H")
text(rep(c(0.0, 0.5), 4), rep(c(1, 0.74, 0.48, 0.2), each = 2), ledders, cex = 1.5)
# dev.off()
# I hope this code and the annotations have been helpfull. If you have any questions please to note hesitate to contact one of us. Beers - CAH & JAF 