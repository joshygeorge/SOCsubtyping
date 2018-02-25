
read.ncounter.normalized.data <- function(fileName)
{

	header <- read.csv(fileName,as.is=T,nrow=15,header=F,row.names=1)
	x <- read.csv(fileName,as.is=T,header=F,skip=15)

	rownames(header) <- gsub(" ", ".", rownames(header))
	rownames(header) <- tolower(rownames(header))
	if ("id" %in% rownames(header))
	{
	        rownames(header)[rownames(header) == "id"] <- "sample.id"
	}

	header <- header[-c(8), ]
	header <- header[, -c(1, 2)]
	header <- header[1:12, ]
	header <- header[, !is.na(header[1, ]) & !is.na(header[2, ])]
	sample.ids <- header["sample.id", ]
	sample.ids <- gsub(" ", ".", sample.ids)
	sample.ids <- gsub("^([0-9])", "X\\1", sample.ids)
	colnames(header) <- sample.ids
	colnames(x) <- c("Code.Class","Name","Accession",sample.ids)
	rownames(x) <- x$Name
	x <- x[,-c(1:3)]
	return(list(exp = x, header = header))
}
