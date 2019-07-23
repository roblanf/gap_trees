library(ggplot2)
library(dplyr)
library(parallel)
library(viridis)

setwd("/Users/roblanfear/Documents/github/gap_trees")
iqpath = "/Users/roblanfear/Documents/github/gap_trees/iqtree"


# set the parameters you like

# internal branch length
int.bl = seq(0.00, 0.010, by=0.001) 

# external branch length
ext.bl = seq(0.00, 0.00, by = 0.01) 

# proportion of gaps
prop_gaps = seq(0.0, 0.99, by = 0.001)

# replicates
reps = seq(1:100)

# force alignments to be variable
forceVarialbe = c(TRUE, FALSE)

# are gaps homologous or not
# homologous is like
# sp_a ------AACGTGCA
# sp_b ------AACGTGCA
homologous = c(TRUE)

# set up a table of all parameters
paramstable = expand.grid("int" = int.bl, "ext" = ext.bl, "gaps" = prop_gaps, "hom" = homologous, "rep" = reps, "forceVarialbe" = forceVarialbe)


# we remove from the paramstable the obviously impossible case
# of monotypic alignments that are forced to be non-monotypic
impossible = which(paramstable$int==0.00 & paramstable$forceVarialbe==TRUE)
paramstable = paramstable[-impossible,]

paramstable$run_ID = 1:nrow(paramstable)

nrow(paramstable)

# get the other functions
source("runIQtree.R")
source("runPhyML.R")
source("run.4tip.gaps.R")
source("pis.R")

# at this point we can consider moving to another folder

# now we get a giant list of the results
tabresults <- mclapply(1:nrow(paramstable), 
					   function(x) run.4tip.gaps(int.bl = paramstable[x, 1], 
												 ext.bl = paramstable[x, 2], 
												 seqlen = 10000, 
												 gaps.prop = paramstable[x, 3], 
												 homgaps = paramstable[x, 4], 
												 savealignment = F, 
												 run_ID = paramstable[x, 7],
												 iqpath = iqpath, 
												 forceNonMonotypic = paramstable[x, 6]),
						mc.cores = 4)

d = as.data.frame(do.call(rbind, tabresults))

d[,c(1:3,6:8,10:11,13:30)] <- apply(d[,c(1:3,6:8,10:11,13:30)], 2, as.numeric)

d$trueTree_search = d$searchTree == "A"

d

save(d, file = "output_data.Rdata")




#### Up to here the analyses seem to work. The plots below have not yet been tested.


e = d %>% group_by(intBl, extBl, gapProp, gapHomol, forceVarialbe) %>% summarise(proportion_true = sum(trueTree_search)/length(trueTree_search))

ggplot(e, aes(x = gapProp, y = proportion_true, colour = as.factor(intBl), group = as.factor(intBl))) + 
	geom_smooth() +
    facet_grid(forceVarialbe~.)

ggplot(e, aes(x = gapProp, y = intBl, fill = proportion_true)) + 
    geom_tile() + 
    scale_fill_viridis() +
    facet_grid(forceVarialbe~.)

ggplot(d, aes(x = intBl, y = parsInfSites)) + geom_smooth() + geom_jitter() + facet_grid(forceVarialbe~.)


