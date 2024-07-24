#Given a vector of statistics and a vector of chromosome labels, 
#this function plots the statistic on separate chromosomes.
#mark.y is a numerical value. This function draws a horiaontal
#line at that value across the whole plot

plot_stats_on_chr <- function(stats, chr, y.label = "Statistic", plot.type = "l", mark.y = NULL){
    
    ymax <- max(stats, na.rm = TRUE)
    ymin <- min(stats, na.rm = TRUE)
    u_chr <- unique(chr)
    n.chr <- length(u_chr)
    
    layout.mat <- matrix(1:(n.chr+2), nrow = 1)
    layout(layout.mat)

    for(p in 1:(n.chr+2)){
        if(p == 1){
            par(mar = c(4,2,4,0))
            plot.new()
            plot.window(xlim = c(0,1), ylim = c(ymin, ymax))
            axis(2, line = -2)
            mtext(y.label, side = 2, line = 0.5)
        }
        if(p == (n.chr+1)){
            par(mar = c(4,0,4,2))
        }
        if(p != 1 && p != (n.chr+2)){
            par(mar = c(4,0,4,0))
            chr.which = p - 1
            chr.idx <- which(chr == chr.which)
            plot.new()
            plot.window(xlim = c(1,length(chr.idx)), ylim = c(ymin, ymax))

            if(p %% 2 == 0){
                draw.rectangle(1, length(chr.idx), ymin, ymax, fill = "gray", border = NA)
            }

            points(stats[chr.idx], type = plot.type)
            mtext(chr.which, side = 1, line = 1)

            if(!is.null(mark.y)){
                abline(h = mark.y)
            }
        }

    }

}