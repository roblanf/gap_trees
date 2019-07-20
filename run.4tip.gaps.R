# This function assigns branch lengths to a 4-tip tree, allowing different internal vs external branch lengths. It simulates an alignment under JC. Two unrelated taxa have a proportion of gaps: (i) at homologous sites, and (ii) at non-homologous sites. It then runs phylogenetic analyses with (i) IQtree, (ii) PhyML, and (ii) MrBayes??, and returns four trees (simulated, and the three estimates).

require(phangorn)

run.4tip.gaps <- function(int.bl = 0.1, ext.bl = 0.1, seqlen = 1000, gaps.prop = 0.0, homgaps = T, savealignment = T, run_ID, iqpath, forceNonMonotypic = F){

	# make tree A
	tr <- read.tree(text = "((t1,t2),t3,t4);")
	tr$edge.length <- c(int.bl, rep(ext.bl, 4))

	# make the alignment
	makeAl <- function(tr, seqlen, gaps.prop, homgaps){
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
	}

	al <- makeAl(tr, seqlen, gaps.prop, homgaps)
	
	while(forceNonMonotypic && length(seg.sites(al)) == 0) al <- makeAl(tr, seqlen, gaps.prop, homgaps)
	write.dna(al, file = paste0(run_ID, ".fasta"), format = "fasta")

	# run IQ-TREE
	iqres <- try(runIQtree(iqpath, run_ID))
	
	if(!savealignment) system(paste0("rm ", run_ID, ".fasta"))

	return(iqres)
}