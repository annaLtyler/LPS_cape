#This function takes in a list of package names.
#It checks to see which are installed already.
#If any are not installed, it installs them, and 
#then loads all packages

load_libraries <- function(package.names, personal.library = FALSE){

    dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)
    
    current.packages <- installed.packages()
    to.install <- setdiff(package.names, rownames(current.packages))

    if(length(to.install) > 0){

        if(personal.library){
            lib.loc <- Sys.getenv("R_LIBS_USER")
        }else{
            lib.loc <- .libPaths()[1]
        }

        for(i in 1:length(to.install)){
            install.packages(to.install[i], lib = lib.loc, 
            repos = "http://cran.us.r-project.org")
        }    

        #check again to see if we missed any. These might be bioconductor packages
        current.packages <- installed.packages()
        to.install <- setdiff(package.names, rownames(current.packages))

        if(length(to.install) > 0){
            if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager", lib = lib.loc, 
            repos = "http://cran.us.r-project.org")

            for(i in 1:length(to.install)){
                BiocManager::install(to.install[i], lib.loc = lib.loc)
            }
        }
    }

    for(i in 1:length(package.names)){library(package.names[i], character.only = TRUE)}

}