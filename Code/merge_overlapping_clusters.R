#sometimes you want clusters to be of adjacent elements
#and you don't want elements clustered with distant 
#clusters. This function does the simplest thing to
#fix this, which is just to merge any clusters that 
#have overlapping elements. 

merge_overlapping_clusters <- function(cl.membership){

    any_overlap <- function(test.cl.idx){
        cl.pairs <- pair.matrix(1:length(test.cl.idx))
        are.overlapping <- rep(NA, nrow = cl.pairs)
        for(clp in 1:nrow(cl.pairs)){
            cl1 <- cl.pairs[clp,1]
            cl2 <- cl.pairs[clp,2]
            are.overlapping[clp] <- as.numeric(segments.overlap(min(test.cl.idx[[cl1]]), 
                max(test.cl.idx[[cl1]]), 
                min(test.cl.idx[[cl2]]), max(test.cl.idx[[cl2]])))
        }
        return(are.overlapping)
    }


    u_cl <- unique(cl.membership)
    cl.idx <- lapply(u_cl, function(x) which(cl.membership == x))
    #plot(cl.membership)
    cl.overlap <- any_overlap(cl.idx)
    merge.cl <- any(which(cl.overlap > 0))

    if(!merge.cl){
        return(cl.membership)
    }

    new.cl <- cl.membership
    while(merge.cl){
        
        u_cl <- unique(new.cl)
        cl.pairs <- pair.matrix(u_cl)

        #just merge any clusters that are overlapping?
        which.overlap <- which(cl.overlap > 0)
        for(clp in which.overlap){
            cl1 <- cl.pairs[clp,1]
            cl2 <- cl.pairs[clp,2]
            cl2.idx <- which(new.cl == cl2)
            if(length(cl2.idx) > 0){
                new.cl[cl2.idx]  <- cl1
            }
        }

        new.u_cl <- unique(new.cl)
        new.cl.labels <- 1:length(new.u_cl)
        renamed.cl <- rep(NA, length(new.cl))
        for(cl in 1:length(new.u_cl)){
            renamed.cl[which(new.cl == new.u_cl[cl])] <- new.cl.labels[cl]
        }

    new.cl <- renamed.cl
    cl.idx <- lapply(new.cl.labels, function(x) which(new.cl == x))
    cl.overlap <- any_overlap(cl.idx)
    merge.cl <- any(which(cl.overlap > 0))
}

    return(new.cl)
    #plot(new.cl)

}