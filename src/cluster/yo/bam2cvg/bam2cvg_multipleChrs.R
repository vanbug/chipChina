# Rscript for finding coverage from BAM files for all chromosomes (Chr1-19,X,Y,M)
# Author : Sukhdeep Singh
# Organization : Max Planck

##########################################

# looping data per chromosome
chip=list();control<-list()

# importing libraries
library(BSgenome.Mmusculus.UCSC.mm9);library(chipseq);library(ShortRead)

# chip and control file user input
#chipFile=scan("stdin",what="character",n=1,quiet=TRUE)
chipFile=commandArgs()[5]
controlFile=commandArgs()[6]

# building chrs list
chrs=c(paste("chr",seq(1,19),sep=''),paste("chr",c('X','Y','M'),sep=''))
what=scanBamWhat()
number=seq(1,22,by=1)

for(i in 1:length(number)){
which<-GRanges((chrs[number[i]]),ranges=IRanges(number[i],c(seqlengths(Mmusculus)[number[i]])))

# replacing NA in GRanges with chr length
seqlengths(which)=seqlengths(Mmusculus)[number[i]]

params<-ScanBamParam(which=which,what=what)

# reading in files
chip=readAligned(chipFile,type="BAM",param=params)
control=readAligned(controlFile,type="BAM",param=params)
print(paste("Read Chr",i,sep=''))

# removing NA's
control=control[which(is.na(control@strand)==F)]
chip=chip[which(is.na(chip@strand)==F)]

# filtering data # FILTERING MADE OFF AS THE PEAK HEIGHT DECREASES SIGNIFICANTLY
#filt<-occurrenceFilter(withSread=FALSE);control<-control[filt(control)];chip<-chip[filt(chip)]

# variable declarations
chip.norm.aln<-list();control.norm.aln<-list();chip.fragment.size<-c();control.fragment.size<-c();chip.cov<-c();control.cov<-c();

# NORMALIZING LIBRARY SIZE FROM CHIP AND CONTROL using sampling
sizeControl<-length(control@id)
sizeChip<-length(chip@id)
smallest.size<-min(sizeControl,sizeChip)

print ("Sampling Datasets")
# normalizing chip and control library size
chip.norm.aln<-sample(chip,smallest.size,replace=FALSE)
control.norm.aln<-sample(control,smallest.size,replace=FALSE)

# fragment length extension using  #SISSR- Site Identification from Short Sequence Reads
chip.fragment.size=round(estimate.mean.fraglen(chip.norm.aln,method="SISSR"))
control.fragment.size=round(estimate.mean.fraglen(control.norm.aln,method="SISSR"))

print ("Obtaining Coverage")

# Coverage
# calculating extend length - length of the read to be extended from 36 bp to 200 or so, to which it was broken down in start
chip.extend.length=chip.fragment.size-width(chip.norm.aln)
control.extend.length=control.fragment.size-width(control.norm.aln)

# calculating coverage, extend - extends the read reference from 36 to 200
chip.cov<-coverage(chip.norm.aln,width=seqlengths(Mmusculus),extend=as.integer(chip.extend.length))
control.cov<-coverage(control.norm.aln,width=seqlengths(Mmusculus),extend=as.integer(control.extend.length))

print (paste("Writing cov for Chr",i))
# writing coverage externally
bedGraph=function(cov,chrs,type){
        pf<-runValue(cov[[1]])>0
        write.table(cbind(chrs,start(cov[[1]])[pf],
        format(end(cov[[1]]),scientific=F)[pf],        
	runValue(cov[[1]])[pf]),file=paste(type,".bedGraph",sep=''),append=TRUE,sep="\t",
        quote=FALSE,row.names=FALSE,col.names=FALSE)
        }

bedGraph(chip.cov,chrs[i],type="chip")
bedGraph(control.cov,chrs[i],type="control")
print (paste("Done with",i))
}
