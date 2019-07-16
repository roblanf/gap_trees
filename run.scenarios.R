run.scenario <- function(row){
	res <- run.4tip.gaps(int.bl = row[1], prop.extbl = row[2], seqlen = 1000, gaps.prop = row[3], homgaps = row[4])
	return(res)      
}

tabresults <- lapply(paramstable, run.scenario)