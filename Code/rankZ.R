#Rank Z normalize values in a vector
#I took this from Dan Gatti's DOQTL package

rankZ <- function (x) 
{
    x = rank(x, na.last = "keep", ties.method = "average")/(sum(!is.na(x)) + 
        1)
    return(qnorm(x))
}
