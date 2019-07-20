# This function extracts some of the main results from an IQtree run output file.

# We return a single-row data frame with the following entries
# 0. The run_ID (this so we can later merge with the parameter table)
# 1. The identity of the tree from the IQ-TREE search ("A", "B", "C", see the treefiles in the repo)
# 2. The lnL of the tree from the IQ-TREE search
# 3. The bootsrap value on the internal branch for the IQ-TREE search
# 4. The tree string for the tree from the IQ-TREE search (in case we want any other information)
# 5. The lnL of the tree_A analysis
# 6. The lnL of the tree_B analysis
# 7. The lnL of the tree_C analysis
# 8-10. The internal branch lengths of the A, B, and C trees
# 11-13. The t1 branch length of the A, B, and C trees
# 14-16. The t2 branch length of the A, B, and C trees
# 17-19. The t3 branch length of the A, B, and C trees
# 20-22. The t4 branch length of the A, B, and C trees 
# 23. The software used to run this, i.e. 'iqtree' (this so we can easily rbind results from different software later)
# 24. The number of variable sites

runIQtree <- function(iqtreePath, run_ID){

    trs <- list(A = read.tree(text = "((t1,t2),t3,t4);"), B = read.tree(text = "((t1,t3),t2,t4);"), C = read.tree(text = "((t1,t4),t3,t2);"))
	
    tr_analysis_ID <- c("", "tree_A", "tree_B", "tree_C")
	
    aln_name = paste0(run_ID, ".fasta")

    # run IQ-TREE without a tree and with 100 bootstraps, then forcing all three resolutions
    
    iqc <- vector()
    iqc[1] <- paste0(iqtreePath, " -s ", aln_name, " -m JC -b 100 ", " -pre ", run_ID)
    iqc[2:4] <- sapply(tr_analysis_ID[2:4], function(x) paste0(iqtreePath, " -s ", aln_name, " -m JC -g ", x, ".tre -pre ", paste0(run_ID, x)))
    for(i in iqc) system(i)
    
    # Extract IQ-TREE output file and tree
    
    iqc_out <- lapply(tr_analysis_ID, function(x) readLines(paste0(run_ID, x, ".iqtree")))
    iqc_trs <- lapply(iqc_out, function(x) read.tree(text = x[grep("Tree in newick format:", x)[1] + 2]))
    
    # Extract parsimony informative sites and number of site patterns

    iqc_parsinf <- as.numeric(gsub(".*:|[(].*| ", "", grep("Number of parsimony informative sites:", iqc_out[[1]], value = T)))
    iqc_sitepat <- as.numeric(gsub(".*:|[(].*| ", "", grep("Number of distinct site patterns:", iqc_out[[1]], value = T)))
    
    # The identity of the tree from the IQ-TREE search ("A", "B", "C", see the treefiles in the repo)
    
    iqc1treeID <- c("A", "B", "C")[sapply(iqc_trs[2:4], function(x) RF.dist(iqc_trs[[1]], x)) == 0]

    # lnL of each of the four analyses
	
    iqc_lik <- sapply(iqc_out, function(x) as.numeric(gsub(".*:|[(].*| ", "", grep("Log-likelihood of the tree:", x, value = T))) )
	
    # Bootstrap support of branch in tree search

    iqc1boot <- iqc_trs[[1]]$node.label[2]
    
    # Tree string from tree search
	
    iqc1trestring <- write.tree(iqc_trs[[1]])
	
    # Branch lengths of A, B, and C trees
	
    iqc_intbl <- sapply(iqc_trs[2:4], function(x) x$edge.length[1])
	
    iqc_extbl <- lapply(iqc_trs[2:4], function(x){
	brlens <- x$edge.length[2:5]
	names(brlens) <- x$tip.label
	brlens <- brlens[c("t1", "t2", "t3", "t4")]
	return(brlens)
    })
    iqc_extbl <- c(sapply(iqc_extbl, function(x) x["t1"]), sapply(iqc_extbl, function(x) x["t2"]), sapply(iqc_extbl, function(x) x["t3"]), sapply(iqc_extbl, function(x) x["t4"]))
    names(iqc_extbl) <- c(paste0(c("A", "B", "C"), "t1len"), paste0(c("A", "B", "C"), "t2len"), paste0(c("A", "B", "C"), "t3len"), paste0(c("A", "B", "C"), "t4len"))
	
    software <- "iqtree"
	
    # Delete unnecessary output
	
    torm <- as.vector(sapply(paste0(run_ID, tr_analysis_ID), function(x) paste0(x, c(".ckp.gz", ".bionj", ".log", ".mldist", ".treefile", ".contree", ".iqtree", ".boottrees", ".parstree"))))
    for(i in torm) system(paste0("rm ", i))
	
    res <- c(runID = run_ID, parsInfSites = iqc_parsinf, sitePats = iqc_sitepat, searchTree = iqc1treeID, searchlnL = iqc_lik[1], searchBoot = iqc1boot, searchTrStr = iqc1trestring, AlnL = iqc_lik[2], BlnL = iqc_lik[3], ClnL = iqc_lik[4], AintBL = iqc_intbl[1], BintBL = iqc_intbl[2], CintBL = iqc_intbl[3], iqc_extbl, software = software)
	
	return(res)

}