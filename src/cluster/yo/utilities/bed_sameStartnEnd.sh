# Rscript which sorts the same start and end postions in a bed file and converts the reference positions of the reads (start and end) of a bed file from scientific 'e' to numerals
# Author : Sukhdeep Singh
# Organization : Max Planck
##############################################################################################################################

bedFile=commandArgs()[5]

# sorts the same start and end read position in a bed file
#covFile=scan("stdin",what="character",n=1,quiet=TRUE)

# reads bed file as table
cov=read.table(bedFile,stringsAsFactors=FALSE)

# boolean matching and adding one the corresponding reference positions
boo=cov[,3]==cov[,2]
cov[,3][boo]<-cov[,3][boo]+1

# we want data in non-scientific format, without 'e'
cov$V2=format(cov$V2,scientific=F)
cov$V3=format(cov$V3,scientific=F)
cov$V4=format(cov$V4,scientific=F)
write.table(cov,file=paste(bedFile,".bedDiff",sep=""),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
print ("Same start and end positions in a bed file sorted. Remember to add +1 in the genome index file while converting bed to bigBed,bigWig or any other format")

#############################################################################################################################
