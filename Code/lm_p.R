#This function extracts a p value from a linear model
lm_p <- function(model){

	f <- summary(model)$fstatistic
    if(is.null(f)){
        return(NA)
        }
	p <- pf(f[1],f[2],f[3],lower.tail=F)
	return(p)
}