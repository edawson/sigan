get.complementary <- function(x) {
	if (x == "A") {
		y <- "T"
	} else if (x == "T") {
		y <- "A"
	} else if (x == "C") {
		y <- "G"
	} else if (x == "G") {
		y <- "C"
	}
	return(y)
}
get.reverse <- function(x) {
	n.x <- length(x)
	z <- rep(NA,4)
	w <- x
	for (i in 1:n.x) {
		y <- as.vector(unlist(strsplit(as.character(x[i]),"")))
		z[1] <- get.complementary(y[1])
		z[2] <- get.complementary(y[2])
		z[3] <- get.complementary(y[4])
		z[4] <- get.complementary(y[3])
		w[i] <- paste(z[1],z[2],z[3],z[4],sep="")
	}
	return(w)
}

context65 <- read.delim(paste(paste("context65.txt",sep=""),sep=""),header=T,sep="\t")
x <- context65[,2]
x1 <- sapply(x,function(x) strsplit(as.character(x)," in ")[[1]][1]) 
x2 <- sapply(x,function(x) strsplit(as.character(x)," in ")[[1]][2])
context65.flanking <- paste(substring(x2,1,1),substring(x2,3,3),sep="")

context96 <- read.delim(paste(paste("context96.txt",sep=""),sep=""),header=F,sep="\t")
context96 <- as.vector(unlist(context96))
context96.reverse <- get.reverse(context96)
context96.label <- c(context96.reverse[1:48],context96[49:96])
context96.label.reverse <- get.reverse(context96.label)
context96.plus <- paste(context96.label,"+",sep="")
context96.minus <- paste(context96.label,"-",sep="")
context96.reverse.plus <- paste(context96.label.reverse,"+",sep="")
context96.reverse.minus <- paste(context96.label.reverse,"-",sep="")

get.indel.spectrum <- function(maf){
    maf[,"sample"] <- maf$Tumor_Sample_Barcode
    if ("Variant_Type"  %in% colnames(maf)){
       maf <- maf[maf$Variant_Type %in% c("INS", "DEL"),]
    }
	if ("SAMPLE_ID" %in% colnames(maf)){
		maf$sample <- maf$SAMPLE_ID
	}
    var.type <- maf$Variant_Type
    indel.base <- toupper(maf[,colnames(maf) %in% c("indel_base")])
    indel.len <- as.numeric(maf[,colnames(maf) %in% c("indel_len")])
    indel.frame <- as.numeric(indel.len) %% 3 == 0
    indel.rpt.count <- maf[,colnames(maf) %in% c("repeat_count")]
    indel.mh.len <- maf[,colnames(maf) %in% c("microhomology_length")]

    contig <- paste(
        var.type,
        indel.base,
        ifelse(indel.len < 4,
		 	indel.len,
			ifelse(indel.len < 10, "4-10", "10+")
		),
		ifelse(indel.len > 3 & indel.mh.len > 1,
		"MH",
		"NoMH")
    , sep = "_")

	# contig <- paste(
	# 	var.type,
	# 	indel.base,
	# 	ifelse(indel.len < 4, as.character(indel.len), "4+"),

	# 	ifelse( indel.len > 1 & (indel.mh.len > 1 | indel.rpt.count > 1),
	# 		ifelse(indel.mh.len > 1, "MH", "RPT"),
	# 		"None"),

	# 	ifelse(indel.len > 1 & (indel.mh.len > 1 | indel.rpt.count > 1),
	# 		ifelse(indel.rpt.count > 1,
	# 			ifelse(indel.rpt.count < 3, indel.rpt.count, "3+"),
	# 			ifelse(indel.mh.len < 3, indel.mh.len, "3+")),
	# 	"X")
	# 	,
	# 	sep = "_"
	# )

	indel.num <- rep(0, length(contig))

	context.indels <- read.delim(paste(paste("ed_indel_context",sep=""),sep=""),header=F,sep="\t")
	context.indels <- as.vector(unlist(context.indels))

	for (i in 1:length(context.indels)){
		indel.num[contig %in% c(context.indels[i])] <- i		
	}
	maf[,"context96.num"] <- NA
	maf[,"context96.word"] <- NA
	maf[,"indel.context"] <- contig
	maf[,"indel.context.num"] <- indel.num

	lego.indels <- as.data.frame.matrix(table(maf$indel.context.num, maf$sample))
	# print(rownames(lego.indels))
	# print(context.indels)
	print(contig[!contig %in% c(context.indels)])

	if (nrow(lego.indels) < length(context.indels)){
		existing <- which(c(1:length(context.indels)) %in% as.numeric(rownames(lego.indels)))
		missing <- which(!c(1:length(context.indels)) %in% as.numeric(rownames(lego.indels)))
		n.missing <- length(missing)
		x.missing <- array(0, dim = c(n.missing, ncol(lego.indels)))
		rownames(x.missing) <- context.indels[missing]
		colnames(x.missing) <- colnames(lego.indels)
		rownames(lego.indels) <- context.indels[existing]
		lego.indels <- rbind(lego.indels, x.missing)
		lego.indels <-lego.indels[match(context.indels, rownames(lego.indels), nomatch = 0),]
	}
	rownames(lego.indels) <- context.indels
	return (list(maf, lego.indels))

}

get.spectrum96.from.maf <- function(maf) {
	maf[,"sample"] <- maf$Tumor_Sample_Barcode
	if ("Variant_Type" %in% colnames(maf)) {
                maf <- maf[maf$Variant_Type %in% c("SNP"),]
    }
	if ("SAMPLE_ID" %in% colnames(maf)){
		maf$sample <- maf$SAMPLE_ID
	}
	ref <- toupper(maf[,colnames(maf) %in% c("Reference_Allele")])
	alt <- toupper(maf[,colnames(maf) %in% c("Tumor_Seq_Allele2")])
	context <- toupper(maf[,colnames(maf) %in% c("ref_context")])
	n.context <- length(unlist(strsplit(as.character(context[1]),"")))
	mid <- trunc(n.context/2)+1
	contig <- paste(ref,alt,substring(context,mid-1,mid-1),substring(context,mid+1,mid+1),sep="")
	context96.num <- rep(0,length(contig))
	for (i in 1:96) {
		context96.num[contig %in% c(context96[i],context96.reverse[i])] <- i
	}
	maf[,"context96.num"] <- context96.num
	maf[,"context96.word"] <- contig
	x <- as.data.frame.matrix(table(maf$context96.num,maf$sample))
	if (nrow(x)<96) {
                existing <- which(c(1:96)%in%as.numeric(rownames(x)))
                missing <- which(!c(1:96)%in%as.numeric(rownames(x)))
                n.missing <- length(missing)
                x.missing <- array(0,dim=c(n.missing,ncol(x)))
                rownames(x.missing) <- context96.label[missing]
                colnames(x.missing) <- colnames(x)
                rownames(x) <- context96.label[existing]
                x <- rbind(x,x.missing)
                x <- x[match(context96.label,rownames(x),nomatch=0),]
	} else if (nrow(x) > 96) {
		stop('Unusal context - double-check context information')
        }
	rownames(x) <- context96.label
	return(list(maf,x))
}


get.spectrum.with.indels <- function(maf){
    #x <- get.spectrum96.from.maf(maf)
	y <- get.indel.spectrum(maf)
	
	## Merge our MAFs and our lego matrices (96 x Nsamples and 128 x Nsamples)
	## First sync up our column names
	#x[[1]][,"indel.context"] <- rep(NA, length(rownames(x)))
	#x[[1]][,"indel.context.num"] <- rep(NA, length(rownames(x)))


	#combined.maf <- rbind(x[[1]], y[[1]])
	#combined.lego <- rbind(x[[2]], y[[2]])
	combined.maf <- y[[1]]
	combined.lego <- y[[2]]
	return (list(combined.maf, combined.lego))

}


