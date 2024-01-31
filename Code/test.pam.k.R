test.pam.k <- function(mat, kseq = 2:10, diss = FALSE, metric = "euclidean", plot.results = TRUE){

    cl.width <- vector(mode = "list", length = length(kseq))
    names(cl.width) <- kseq

    mem <- matrix(NA, nrow = nrow(mat), ncol = length(kseq))
    colnames(mem) <- kseq
    rownames(mem)  <- rownames(mat)

    for(i in 1:length(kseq)){
        km <- pam(mat, k = kseq[i], diss = FALSE, metric = metric)
        mem[,i] <- km$clustering
        mem.f <- data.frame(as.factor(mem[,i]))
        si <- silhouette(km)
        if(plot.results){plot(si)}
        cl.width[[i]] <- summary(si)$clus.avg.widths
    }
    return(list("cl.width" = cl.width, "mem" = mem))
}