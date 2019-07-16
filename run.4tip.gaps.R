# This function assigns branch lengths to a 4-tip tree, allowing different internal vs external branch lengths. It simulates an alignment under JC. Two unrelated taxa have a proportion of gaps: (i) at homologous sites, and (ii) at non-homologous sites. It then runs phylogenetic analyses with (i) IQtree, (ii) PhyML, and (ii) MrBayes??, and returns four trees (simulated, and the three estimates).

require(phangorn)

run.4tip.gaps <- function(int.bl = 0.01, prop.extbl = 1, seqlen = 1000, gaps.prop = 0.1, homgaps = T, tempname = "temp.alignment.fasta", iqpath = "/Users/roblanfear/Dropbox/Projects_Current/magic_trees/iqtree", phymlpath = "/Users/roblanfear/Dropbox/Projects_Current/magic_trees/iqtree", savealignment = T){
	tr <- read.tree(text = "((t1,t2),t3,t4);")
	tr$edge.length <- c(int.bl, rep(int.bl*prop.extbl, 4))
	al <- as.character(as.matrix(as.DNAbin(simSeq(tr, l = seqlen))))
	if(gaps.prop > 0){
		ngapsites <- round(seqlen * gaps.prop)
		if(homgaps){
			al[c("t1", "t3"), 1:ngapsites] <- "-"
		} else {
		        al["t1", 1:ngapsites] <- "-"
			al["t3", ncol(al):(ncol(al)-ngapsites+1)] <- "-"
		}
	}
	al <- as.DNAbin(al)
	write.dna(al, file = tempname, format = "fasta")
	#return(al)
	temptree <- tr
	temptree$edge.length <- NULL

	iqres <- try(runIQtree(tempname, format = 'fasta', aadata = F, tempname, iqtreePath = iqpath, model = 'JC'))
	
	#phymlres <- try(runPhyML(tempname, format = 'fasta', aadata = F, tempname, phymlPath = phymlpath, model = 'JC'))

	trees = list(tr, iqres$tree)
	class(trees) <- "multiPhylo"

	print(trees[[1]])
	print(trees[[2]])

	if(!savealignment) system(paste0("rm ", tempname))
	res <- list(simulated.tree = tr, iqtree.results = iqres, topdist = RF.dist(trees))
	return(res)
}