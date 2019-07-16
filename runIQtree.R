# This function extracts some of the main results from an IQtree run output file.

runIQtree <- function(iqtreePath, run_ID){

    aln_name = paste0(run_ID, ".fasta")

    # run IQtree without a tree and with 100 bootstraps
    iqc1 = paste0(iqtreePath, " -s ", aln_name, " -m JC -alrt 1000 -b 100 ", " -pre ", run_ID)
    system(iqc1)

    # run IQtree forcing all three resolutions
    iqc2 = paste0(iqtreePath, " -s ", aln_name, " -m JC -g tree_A", treefile, " -pre ", paste0(run_ID, "tree_A"))
    system(iqc2)
    iqc3 = paste0(iqtreePath, " -s ", aln_name, " -m JC -g tree_B", treefile, " -pre ", paste0(run_ID, "tree_B"))
    system(iqc3)
    iqc4 = paste0(iqtreePath, " -s ", aln_name, " -m JC -g tree_C", treefile, " -pre ", paste0(run_ID, "tree_C"))
    system(iqc4)
    
    
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





    # get the tree from the search
    allout <- readLines(paste0(fileName, ".iqtree"))

  	torm <- grep(paste0(fileName, c(".ckp.gz", ".bionj", ".log", ".mldist", ".treefile", ".uniqueseq.", ".contree", ".splits.nex", ".iqtree", ".parstree"), collapse = "|"), dir(), value = T)
  	
  	for(i in torm) system(paste0("rm ", i))
         
      res <- list()

      esttree <- read.tree(text = allout[grep("Tree in newick format:", allout)[1] + 2])

      boots <- as.numeric(esttree$node.label)
      res$meanNodeSupport <- mean(boots, na.rm = T)

      res$nodeSupport95 <- diff(quantile(boots, c(0.025, 0.975), na.rm = T))

      trlen <- grep("Total tree length", allout, value = T)
      res$treeLength <- as.numeric(gsub(".* ", "", trlen))

      if(length(grep("HKY|GTR", model))==1){
      	piLocation <- grep('State frequencies', allout)+2
      	piParams <- as.numeric(gsub('.* ', '', allout[(piLocation):(piLocation+3)]))
      	res$piParams <-  piParams
      }

      if(length(grep("[+]G", model))==1){
      	alphaLocation <- grep('Gamma shape alpha', allout)
      	alphaParam <- as.numeric(gsub(".* ", "", allout[alphaLocation]))
      		res$alphaParam <- alphaParam
      }

      if(length(grep("GTR", model))==1){
      	gtrLocation <- grep('Rate parameter R', allout)+2
      	gtrMatrix <- as.numeric(gsub('.* ', '',  allout[(gtrLocation):(gtrLocation+5)]))
      	res$gtrMatrix <-  gtrMatrix
      }

      if(length(grep("HKY", model))==1){
      	trtvLocation <- grep('A-G:', allout)
      	trtvRatio <- as.numeric(gsub(".* ", "", allout[trtvLocation]))
      	res$trtvRatio <-  trtvRatio
      }

      likse <- grep("Log-likelihood of the tree:", allout, value = T)
      res$likelihood <- as.numeric(gsub(".*:|[(].*| ", "", likse))

      unclik <- grep("Unconstrained log-likelihood", allout, value = T) 
      res$uncLikelihood <- as.numeric(gsub(".* ", "", unclik))

      res$delta <- res$uncLikelihood - res$likelihood

      res$tree <- esttree

      # Other parameters that might be of interest in the future
      #params <- grep("Number of free parameters", allout, value  = T)
      #res$Nparams <- gsub(".* ", "", params)       
      #aic <- grep("Akaike information criterion", allout, value = T)
      #res$aic <- gsub(".* ", "", aic[2])
      #res$aicc <- gsub(".* ", "", aic[3])       
      #bic <- grep("Bayesian information criterion", allout, value = T)
      #res$bic <- gsub(".* ", "", bic[2])       
      #bootrange <- range(boots, na.rm = T)
      #res$bootmax <- max(bootrange)
      #res$bootmin <- min(bootrange)
      #res$bootrangesize <- diff(bootrange)       
      #res$likelihoodSE <- gsub(".*[(]|[)].*|.* ", "", likse)

      return(res)
}