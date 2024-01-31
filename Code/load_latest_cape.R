#This function loads a local version of cape code.

load_latest_cape <- function(cape.dir, personal.library = FALSE){

    needed.libraries <- c("here", "qtl2", "abind", "Matrix", "MASS", "regress", "igraph",
        "RColorBrewer", "R6", "yaml", "tools", "caTools", "propagate", "igraph",
        "qtl", "regress", "evd", "shape", "pracma", "qtl2convert")
    load_libraries(needed.libraries, personal.library = personal.library)

    cape.fun <- list.files(file.path(cape.dir, "R"), pattern = ".R", full.names = TRUE)
    for(i in 1:length(cape.fun)){source(cape.fun[i])}

}