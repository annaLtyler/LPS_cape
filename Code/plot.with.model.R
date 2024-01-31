#This function takes an x and y and plots them with a fitted linear model
#if write.results is FALSE, the correlation or model results will not be
#written in the plot title

plot.with.model <- function(x, y, xlim = NULL, ylim = NULL, col = "black", 
line.col = "#a6bddb", pch = 16, main, xlab, ylab, report = c("lm", "cor.test"), 
cex = 1, plot.results = TRUE, write.results = TRUE, add = FALSE, plot.type = "p"){
	
	if(missing(main)){main = ""}
	if(missing(xlab)){xlab = deparse(substitute(x))}
	if(missing(ylab)){ylab = deparse(substitute(y))}

	report <- report[1]
	
	model <- lm(y~x)

	if(!is.na(coef(model)[2])){
		if(report == "lm"){
			f <- summary(model)$fstatistic
			p <- pf(f[1],f[2],f[3],lower.tail=F)
			r2 <- signif(summary(model)$r.squared, 2)
			rlab <- "R2"
			}else{
			test <- cor.test(x,y)	
			p <- signif(test$p.value, 2)
			r2 <- signif(test$estimate, 2)
			rlab <- "r"
			}
	}else{
		if(report == "lm"){
			f <- NA
			p <- NA
			r2 <- NA
			rlab <- "R2"
		}else{
			test <- NA
			p <- NA
			r2 <- NA
			rlab <- "r"
			}
	}

	if(plot.results){
		if(write.results){
			new.title <- paste0(main, "\n", rlab, " = ", r2, ", p = ", signif(p, 2))
		}else{
			new.title <- ""
		}
	if(!add){
		plot(x,y, xlim = xlim, ylim = ylim, col = col, pch = pch, main = new.title, 
		xlab = xlab, ylab = ylab, type = plot.type, cex = cex)
	}else{
		points(x,y, col = col, pch = pch, cex = cex)
	}

	if(!is.na(coef(model)[2])){
		abline(model, col = line.col, lwd = 3)
	}
	}
	
	invisible(c("r" = r2, "p" = p))
	
}