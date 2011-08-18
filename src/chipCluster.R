###################################################################################################################################
# R script for analysis of Chip Seq data after Shell pre-processingaligned data
# Author : Sukhdeep Singh
# Organization : Max Planck Dresden
###################################################################################################################################

# importing libraries
library('ShortRead');
library('chipseq');

# declaring variables
chrs<-c();chipChr<-list();controlChr<-list();

# reading the user inputted file
control<-readline('Please enter the control location! ')
control<-gsub('\'','',control)
chip<-readline('Please enter the chip location! ')
chip<-gsub('\'','',chip)


# reading control and chip data
control<-readAligned(control,type="BAM");
chip<-readAligned(chip,type="BAM");

# printing info
print (control);
print (chip);

# uniquely aligned data, dropping NA's
controlP=control[which(is.na(control@strand)==F)]
chipP=chip[which(is.na(chip@strand)==F)]

# inbuilt filtering for alignQuality occurrenceFilter
filt1<-chromosomeFilter("chr[0-9XYM]")
filt2<-alignQualityFilter(25)
filt3<-occurrenceFilter(withSread=FALSE)
filt<-compose(filt1,filt2,filt3)
controlF<-controlP[filt(controlP)]
chipF<-chipP[filt(chipP)]

# sort chromosome levels 
sortChr=function(x){
inputlevel=levels(x)
sortedchrlevel=paste("chr",sort(as.numeric(gsub("chr","",as.character(inputlevel)))),sep="")
diff=length(inputlevel)-length(sortedchrlevel)
if (length(inputlevel)>length(sortedchrlevel)&&(diff!=3)) {print ("Investigate")} 
else {sortedchrlevel[(length(sortedchrlevel)+1)]="chrX";sortedchrlevel[(length(sortedchrlevel)+1)]="chrY";sortedchrlevel[(length(sortedchrlevel)+1)]="chrM"}
sortedchrlevel=paste(sortedchrlevel,"$",sep='')
#sortedchrlevel[(length(sortedchrlevel)+1):(length(sortedchrlevel)+3)]=c("chrX","chrY","chrM")
return (sortedchrlevel)
}

# simplified chromosome extraction - single chromosome, multiple sample extraction
extractChr=function(control1level=NULL,control2level=NULL,chip1level=NULL,chip2level=NULL,chip3level=NULL,control1=NULL,control2=NULL,chip1=NULL,chip2=NULL,chip3=NULL){
if (is.null(control1)==FALSE) {controlChr1=control1[(chromosomeFilter(control1level))(control1)]} else {controlChr1=NULL}
if (is.null(control2)==FALSE) {controlChr2=control2[(chromosomeFilter(control2level))(control2)]} else {controlChr2=NULL}
if (is.null(chip1)==FALSE)    {chipChr1=chip1[(chromosomeFilter(chip1level))(chip1)]} else {chipChr1=NULL}
if (is.null(chip2)==FALSE) {chipChr2=chip2[(chromosomeFilter(chip2level))(chip2)]} else {chipChr2=NULL}
if (is.null(chip3)==FALSE) {chipChr3=chip3[(chromosomeFilter(chip3level))(chip3)]} else {chipChr3=NULL}
return (list(control1=controlChr1,control2=controlChr2,chip1=chipChr1,chip2=chipChr2,chip3=chipChr3))
}

# lappy usage - multiple chromosome, multiple sample extraction
extractNucleosome=function(control1level=NULL,chip1level=NULL,control1=NULL,chip1=NULL){
chrs<-mapply(extractChr,control1level=control1level,chip1level=chip1level,MoreArgs=list(control1=control1,chip1=chip1))
return (chrs)
}

# sort chromosome levels
#controlLevels=sortChr(control@chromosome)
#chipLevels=sortChr(chip@chromosome)

# extract chip's n control's chromosome specific data
#extractedChrs<-extractNucleosome(controlLevels,chipLevels,controlP,chipP)


#fetching chromosome data using seq
#controlChrs=extractedChrs[seq(1,length(extractedChrs),by=5)]
#chipChrs=extractedChrs[seq(3,length(extractedChrs),by=5)]

# unlisting chip and control chromosome data to remove head list - NO NEED WITH CURRENT CODE
#controlChr=unlist(controlChr)
#chipChr=unlist(chipChr)

# variable declarations
chip.norm.aln<-list(); control.norm.aln<-list()
chip.fragment.size<-c(); control.fragment.size<-c()
chip.cov<-c(); control.cov<-c()

# NORMALIZING LIBRARY SIZE FROM CHIP AND CONTROL - in GFP & MEN, control has more reads
for (i in 1:length(chipChrs)){
sizeControl<-length(controlChrs[[i]]@id)
sizeChip<-length(chipChrs[[i]]@id)
smallest.size<-min(sizeControl,sizeChip)
chip.norm.aln[[i]]<-sample(chipChrs[[i]],smallest.size,replace=FALSE)
control.norm.aln[[i]]<-sample(controlChrs[[i]],smallest.size,replace=FALSE)
print (paste("Done with Chr",i,sep=''))
#print (chip.norm.aln[[i]])
#print (control.norm.aln[[i]])
chip.fragment.size=c(chip.fragment.size,round(estimate.mean.fraglen(chip.norm.aln[[i]],method="SISSR")))
control.fragment.size=c(control.fragment.size,round(estimate.mean.fraglen(control.norm.aln[[i]],method="SISSR")))
}

# Investigating Coverage
# fetching the sequence lengths of mm9 - the chromosome size
library(BSgenome.Mmusculus.UCSC.mm9)
seqlengths(Mmusculus)

# calculating extend length
chip.extend.length=chip.fragment.size-width(chip.norm.aln[[1]][1])
control.extend.length=control.fragment.size-width(control.norm.aln[[1]][1])

# calculating coverage
for (i in 1:length(chipChrs)){
chip.cov<-c(chip.cov,coverage(chip.norm.aln[[i]],width=seqlengths(Mmusculus),extend=as.integer(chip.extend.length[i])))
control.cov<-c(control.cov,coverage(control.norm.aln[[i]],width=seqlengths(Mmusculus),extend=as.integer(control.extend.length[i])))
print (paste((length(chipChrs)-i),"chromosomes left",sep=''))
}

# JUST WITH INITIAL ANALYSIS - NO NEED OF RECURRANCE
# plotting function for coverage visualization of each chromosome
plotChIP.Coverage<-function(x,xlab="Position",ylab="Coverage",main="ChIP Coverage",sub){
	plot(c(start(x),length(x)),c(runValue(x),1),type="l",col="blue",xlab=xlab,ylab=ylab,main=main,sub=sub)
}

plotControl.Coverage<-function(x,xlab="Position",ylab="Coverage",
		main="Control Coverage",sub)
{
	plot(c(start(x),length(x)),c(runValue(x),1),type="l",
			col="red",xlab=xlab,ylab=ylab,main=main,sub=sub)
}

# plots coverage for control and chip single or multiple chromosomes, x=chip.cov
plotChr=function(x,type){
	for (i in 1:length(x)) {
	if (type=="control") {jpeg(paste("chr",i,".jpg",sep=''))
				plotControl.Coverage(x[[i]][[1]],sub=paste(type,"chr",i))} else {
				plotChIP.Coverage(x[[i]][[1]])}
	dev.off()
	}
}

# coverage to track data
tracks=function(x){return (as(x,"RangedData"))}
nozero=function(x){return (subset(x,score>0))}
# obtaining tracks
chip.cov.Track=lapply(chip.cov,tracks)
control.cov.Track=lapply(control.cov,tracks)
chip.cov.no.zero=lapply(chip.cov.Track,nozero)
control.cov.no.zero=lapply(control.cov.Track,nozero)

# exporting the coverage no zero files as bedGraph files for each chromosome
exportBed=function(cov.no.zero,bedFileName){
export(cov.no.zero,bedFileName,"bedGraph")
}
# have to run a loop instead of mapply so as to get distinct names for different chrmosome bedFiles
for (i in 1:length(control.cov.no.zero)){
exportBed(control.cov.no.zero[[i]],paste(controlLevels[i],".control",sep=''))
exportBed(chip.cov.no.zero[[i]],paste(chipLevels[i],".chip",sep=''))
}

# single line production of bed files is, distinct filename feature is sorted.
#mapply(exportBed,chip.cov.no.zero,list=(bedFileName="chip"))
##########################################
# End of Code

