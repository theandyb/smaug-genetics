##############################################################################
# Add Bin/Category/Sequence columns to summary file, update colnames of bin
# file, and merge overall count for each sequence motif into summary file
##############################################################################
updateData <- function(summfile, binfile, adj){

	nbp <- adj*2+1

	# summfile$BIN <- ceiling(summfile$POS/binw)
	cat("Assigning summfile categories...\n")
	summfile$CAT <- paste(summfile$REF, summfile$ALT, sep="")

	# Manually remove bins near chr20 centromere
	# chr22 <- chr22[ which(chr22$BIN<260 | chr22$BIN>300),]
	summfile$Category[summfile$CAT=="AC" | summfile$CAT=="TG"] <- "AT_CG"
	summfile$Category[summfile$CAT=="AG" | summfile$CAT=="TC"] <- "AT_GC"
	summfile$Category[summfile$CAT=="AT" | summfile$CAT=="TA"] <- "AT_TA"
	summfile$Category[summfile$CAT=="GA" | summfile$CAT=="CT"] <- "GC_AT"
	summfile$Category[summfile$CAT=="GC" | summfile$CAT=="CG"] <- "GC_CG"
	summfile$Category[summfile$CAT=="GT" | summfile$CAT=="CA"] <- "GC_TA"

	cat("Assigning summfile sequences...\n")
	summfile$Sequence <- ifelse(
		substr(summfile$SEQ,adj+1,adj+1)<substr(summfile$ALTSEQ,adj+1,adj+1),
		paste0(summfile$SEQ,"(",summfile$ALTSEQ,")"),
		paste0(summfile$ALTSEQ,"(",summfile$SEQ,")")
	)

	summfile$SEQMIN <- pmin(summfile$SEQ, summfile$ALTSEQ)

	cat("Assigning summfile CpG categories...\n")
	# Second category column to include +3 CpG categories
	summfile$Category2 <- ifelse(substr(summfile$Sequence,adj+1,adj+2)=="CG",
								paste0("cpg_",summfile$Category),
								summfile$Category)

	# get complement of sequence columns in bin file and remove duplicates
	# change this to an apply()-based function?
	cat("Updating bin file columns...\n")
	for(i in 6:ncol(binfile)){
		names(binfile)[i] <- paste0(names(binfile)[i], "(", revcomp(names(binfile)[i]), ")" )
	}

	cat("Removing redundant columns from bin file...\n")
	bins2 <- binfile[,names(binfile)%in%unique(summfile$Sequence)]
	binfile <- cbind(binfile[,1:5],bins2)
	xmax <- floor(max(summfile$BIN)/100)*100

	cat("Counting total motifs in genome...\n")
	mct <- melt(binfile[,5:ncol(binfile)], id="BIN") %>%
		group_by(variable) %>%
		summarise(value=sum(value))
	# bins2 <- aggregate(data=bins2, value ~ variable, sum)
	names(mct) <- c("Sequence", "COUNT")
	mct$Sequence <- sub("[.]", "(", mct$Sequence)
	mct$Sequence <- sub("[.]", ")", mct$Sequence)
	mct$SEQ1 <- substr(mct$Sequence, 0, adj*2+1)
	mct$SEQ2 <- substr(mct$Sequence, (adj*2+1)+2, (adj*2+2)+(adj*2+1))
	mct$SEQMIN <- pmin(mct$SEQ1, mct$SEQ2)
	mct <- mct %>%
		dplyr::select(SEQMIN, COUNT)
	# Find a more efficient way to do this?
	# cat("Updating summary file with motif counts...\n")

	# summfile <- left_join(summfile, bins2, by="SEQMIN")

	datalist<- list("summ"=summfile, "bin"=binfile, "mct"=mct)
	return(datalist)
}