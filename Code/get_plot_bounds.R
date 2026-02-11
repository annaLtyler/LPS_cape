#This function helps decide reasonable plotting boundaries if you
#have a fractional min/max, and you want to extend the boundary
#just beyound those values. The function multiplies your values
#by a factor (scale.factor), rounds up or down as appropriate,
#and then divides by the same factor.
#if center is TRUE, the values will be centered around 0.
#if you want to make sure the returned values end in an odd or an 
#even digit, set return.even or return.odd to TRUE. The default is
#to return a number ending in an even digit.


get_plot_bounds <- function(min.val, max.val, scale.factor = 100, return.even = TRUE, 
    return.odd = FALSE, center = FALSE){

    if(return.even && return.odd){stop("Only one of return.even or return.odd can be set to TRUE.")}

    rounded.min <- floor(min.val*scale.factor)
    rounded.max = ceiling(max.val*scale.factor)

    if(return.even){
        if(abs(rounded.min) %% 2 == 1){
        rounded.min = rounded.min - 1
        }
        if(abs(rounded.max) %% 2 == 1){
        rounded.max = rounded.max + 1
        }
    }

    if(return.odd){
        if(abs(rounded.min) %% 2 == 0){
        rounded.min = rounded.min - 1
        }
        if(abs(rounded.max) %% 2 == 1){
        rounded.max = rounded.max + 1
        }
    }

    if(center){
        max.abs <- max(abs(c(rounded.min, rounded.max)))
        final.min <- (max.abs/scale.factor)*-1
        final.max <- (max.abs/scale.factor)
    }else{
        final.min <- rounded.min/scale.factor
        final.max <- rounded.max/scale.factor
    }

    return(c(final.min, final.max))
}