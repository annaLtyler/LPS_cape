
#VCF processing from copilot

#install.packages("vcfR")
#install.packages("dplyr")

setwd("~/Documents/Projects/LPS/LPS_cape")

library(here)
library(vcfR)
library(dplyr)
library(gprofiler2)
vcf_files <- list.files(here("Data", "vcf"), full.names = TRUE)


#read_vcf <- function(file) {
#  vcf <- read.vcfR(file)
#  vcf_df <- as.data.frame(vcf@fix)
#  return(vcf_df)
#}

extract_region <- function(vcf_file, chrom, start.pos, end.pos) {
  vcf <- read.vcfR(vcf_file)
  chrom.idx <- which(vcf@fix[,"CHROM"] == chrom)
  after.start <- which(vcf@fix[,"POS"] >= start.pos)
  before.end <- which(vcf@fix[,"POS"] <= end.pos)
  region.idx <- Reduce("intersect", list(chrom.idx, after.start, before.end))
  region <- as.data.frame(vcf@fix[region.idx,])
  return(region)
}

#Fshr region: chr17:88,617,268-89,109,594
chr = 17; start.pos = 88617268; end.pos = 89109594 #region around Lhcgr and Fshr
chr = 17; start.pos = 86946507; end.pos = 89246304; #larger region including Cript, an RNA splicing gene

region <- lapply(vcf_files, function(x) extract_region(x, chr, start.pos = start.pos, end.pos = end.pos))

table.names <- gsub(".mgp.v5.indels.dbSNP142.normed.vcf.gz", "", basename(vcf_files))

merged_region <- region[[1]]
for(i in 2:length(region)){
	merged_region <- merge(merged_region, region[[i]], by = c("CHROM", "POS"), 
		suffixes = c("", paste0("_", table.names[i])))
}
colnames(merged_region)[3:ncol(region[[1]])] <- paste(colnames(merged_region)[3:ncol(region[[1]])], table.names[1], sep = "_")


base.col <- colnames(region[[1]])[3:ncol(region[[1]])]
col.idx <- lapply(base.col, function(x) grep(x, colnames(merged_region)))
names(col.idx) <- base.col

alt.table <- ref.table <- matrix(NA, nrow = nrow(merged_region), ncol = length(table.names))
colnames(alt.table) <- colnames(ref.table) <- table.names
rownames(alt.table) <- merged_region[,"POS"]
n.pos <- nrow(merged_region)
for(p in 1:n.pos){
	chr.pos <- merged_region[p,"POS"]
	ref <- as.character(merged_region[p,col.idx$REF])
	alt <- as.character(merged_region[p,col.idx$ALT])
	#the reference and alternate alleles seem to be mixed up.
	#adjust so the reference allele is the shorter of the two strings
	test.table <- cbind(ref, alt)
	min.len <- t(apply(test.table, 1, function(x) which.min(nchar(x))))
	max.len <- t(apply(test.table, 1, function(x) which.max(nchar(x))))
	new.ref <- sapply(1:length(min.len), function(x) test.table[x,min.len[x]])
	new.alt <- sapply(1:length(max.len), function(x) test.table[x,max.len[x]])
	ref.table[p,] <- new.ref
	alt.table[p,] <- new.alt
}

#find out which positions have variation
pos.var <- which(apply(alt.table, 1, function(x) !all(x == x[1])))
write.table(alt.table[pos.var,], here("Results", "VCF", "Variant_Table.tsv"), quote = FALSE,
	sep = "\t")

info.table <- as.matrix(merged_region[pos.var, col.idx$INFO])

make_info_table <- function(info.row){
	split.info <- strsplit(info.row, "|", fixed = TRUE)
	num.col <- sapply(split.info, length)
	split.table <- matrix(NA, nrow = length(table.names), ncol = max(num.col))
	for(i in 1:length(split.info)){
		split.table[i,1:num.col[i]] <- split.info[[i]]
	}
	rownames(split.table) <- table.names
	return(split.table)
}

get_csq <- function(var.table){
	csq.sep <- sapply(strsplit(var.table[,1], ";"), function(x) x[4])
	int.csq <- which(csq.sep != "CSQ=-")
	if(length(int.csq) == 0){
		return(NULL)
	}
	
	strain.csq <- vector(mode = "list", length = length(int.csq))
	names(strain.csq) <- paste(names(int.csq), csq.sep[int.csq], sep = "_")
	
	for(i in 1:length(int.csq)){
		gene.idx <- grep("ENSMUSG", var.table[int.csq[i],])
		if(length(gene.idx) > 0){
			transcript.idx <- grep("ENSMUST", var.table[int.csq[i],])
			gene.id <- var.table[int.csq[i],gene.idx]
			gene.name <- gconvert(gene.id, organism = "mmusculus")[,"name"]
			trans.id <- var.table[int.csq[i],transcript.idx]
			all.csq <- var.table[int.csq[i],(transcript.idx+2)]
			csq.table <- cbind(gene.id, gene.name, trans.id, all.csq)
			strain.csq[[i]] <- csq.table
			}else{
				intergenic.idx <- which(var.table[int.csq[i],] == "intergenic_variant")
				if(length(intergenic.idx) > 0){
					strain.csq[[i]] <- "intergenic_variant"
					}
			}
		}

	return(strain.csq)
}


test <- apply(info.table, 1, make_info_table)
csq <- lapply(test, get_csq)
names(csq) <- rownames(alt.table)[pos.var]

has.vals <- which(sapply(csq, length) > 0)
has.csq <- csq[has.vals]

info.name <- here("Results", "VCF", "Variant_Info.tsv")
sink(info.name)
for(i in 1:length(has.csq)){
	cat("Position:", names(has.csq)[i], "\n")
	for(j in 1:length(has.csq[[i]])){
		cat(names(has.csq[[i]])[j], "\n")
		if(class(has.csq[[i]][[j]])[1] == "character"){
				cat(has.csq[[i]][[j]], "\n")
			}else{
				apply(has.csq[[i]][[j]], 1, function(x) cat(x, "\n"))
			}
	}
	cat("\n\n")
}
sink()

csq.genes <- NULL
for(i in 1:length(has.csq)){
	if(length(has.csq[[i]]) > 0){
		for(j in 1:length(has.csq[[i]])){
			if(length(has.csq[[i]][[1]]) > 1){
				csq.genes <- c(csq.genes, has.csq[[i]][[j]][,"gene.name"])
			}
		}
	}
}

sort(unique(csq.genes))

#========================================================================
#========================================================================
# Read VCF files for each strain
#vcf_data <- lapply(vcf_files, read_vcf)

# Merge dataframes on CHROM and POS
merged_vcf <- Reduce(function(x, y) merge(x, y, by = c("CHROM", "POS"), suffixes = c("", "_other")), vcf_data)

# Save the merged dataframe to a new CSV file
write.csv(merged_vcf, 'merged_vcf.csv', row.names = FALSE)