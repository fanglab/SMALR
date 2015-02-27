require("seqinr")
library(tools)

transMotif<-function(motif)	{
	motif<-gsub("W","[AT]",motif)
	motif<-gsub("S","[CG]",motif)
	motif<-gsub("M","[AC]",motif)
	motif<-gsub("K","[GT]",motif)
	motif<-gsub("R","[AG]",motif)
	motif<-gsub("Y","[CT]",motif)
	motif<-gsub("B","[CGT]",motif)
	motif<-gsub("D","[AGT]",motif)
	motif<-gsub("H","[ACT]",motif)
	motif<-gsub("V","[ACG]",motif)
	motif<-gsub("[N-]","[ACGT]",motif)
	
	motif<-gsub("A","a",motif)
	motif<-gsub("T","t",motif)
	motif<-gsub("G","g",motif)
	motif<-gsub("C","c",motif)
	
	return(motif)
}

strReverse <- function(x)
        sapply(lapply(strsplit(x, NULL), rev), paste, collapse = "")

reverseMotif<-function(motif) {
	motif<-gsub("a","1",motif)
	motif<-gsub("g","2",motif)
	motif<-gsub("c","3",motif)
	motif<-gsub("t","4",motif)
	motif<-gsub("\\[","5",motif)
	motif<-gsub("\\]","6",motif)
	
	motif<-gsub("1","t",motif)
	motif<-gsub("2","c",motif)
	motif<-gsub("3","g",motif)
	motif<-gsub("4","a",motif)
	motif<-gsub("5","\\]",motif)
	motif<-gsub("6","\\[",motif)
	
	motif=strReverse(motif)
	return(motif)
}


findMotifReadSpace<-function(motif,shift,genomeName) {
	
	motifF<-transMotif(motif)
	motifB<-reverseMotif(motifF)
	motifF_Pos<-shift
	motifB_Pos<-nchar(motif)-shift-1
	label<-file_path_sans_ext(basename(genomeName))

	genome<-read.fasta(file = genomeName, as.string = TRUE, seqtype = "DNA")[[1]]	
	
		position_pos<-gregexpr(motifF,genome)[[1]]+motifF_Pos-1
		position_pos_name<-paste(label, "_", motif, "_pos.txt",sep="")
		print(position_pos_name)
		write.table(position_pos,position_pos_name,col.name=FALSE,row.name=FALSE)
		
		position_neg<-gregexpr(motifB,genome)[[1]]+motifB_Pos-1
		position_neg_name<-paste(label, "_", motif, "_neg.txt",sep="")
		print(position_neg_name)
		write.table(position_neg,position_neg_name,col.name=FALSE,row.name=FALSE)
		

	
}

