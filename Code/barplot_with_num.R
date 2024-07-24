#This barplotting function plots numbers on
#top of each bar.


barplot_with_num <- function(values, xlab = "", ylab = "", main = "", col = "gray",
    cex = 1, text.gap = 0.1, text.shift = 0.05, text.srt = 0, horiz = FALSE, 
    names = NULL, las = 1, ylim = NULL, at = NULL, adj = 0.5){

    plot.range <- max(values, na.rm = TRUE) - min(values, na.rm = TRUE)
    coords <- barplot(values, plot = FALSE)

    text.pos <- values+(plot.range*text.shift)
    neg.idx <- which(values < 0)
    if(length(neg.idx) > 0){
        text.pos[neg.idx] <- values[neg.idx]-(plot.range*text.shift)
    }

    if(!horiz){
        if(is.null(ylim)){
            ylim <- c(min(values)-(plot.range*text.gap), max(values)+(plot.range*text.gap))
            if(ylim[1] > 0){ylim[1] <- 0}
        }
        xlim <- c(0, (max(coords)+0.5))

        a <- barplot(values, xlab = xlab, ylab = ylab, main = main, col = col, 
        ylim = ylim, xlim = xlim, horiz = horiz, names = names, las = las)
        text(x = a[,1], y = text.pos, labels = values, cex = cex, srt = text.srt, adj = adj)

    }
    if(horiz){
        if(is.null(ylim)){
            xlim <- c(min(values)-(plot.range*text.gap), max(values)+(plot.range*text.gap))
        }
        ylim <- c(0, (max(coords)+0.5))   

        a <- barplot(values, xlab = xlab, ylab = ylab, main = main, col = col, 
        ylim = ylim, xlim = xlim, horiz = horiz, names = names, las = las)
        text(x = text.pos, y = a[,1], labels = values, cex = cex, at = at, adj = adj)
    }

}