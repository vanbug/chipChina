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
sortedchrlevel=paste(l,"$",sep='')
#sortedchrlevel[(length(sortedchrlevel)+1):(length(sortedchrlevel)+3)]=c("chrX","chrY","chrM")
return (sortedchrlevel)
}

# chromosome specific control data filtering using ShortRead filter - strange bug for chr-1 (picks all with 1 and doublets-10,11 etc) #EDIT : use $ to end the grep match
extChr_Control=function(x){
chrs=paste("chr",x,"$",sep="")
y=controlF[(chromosomeFilter(chrs))(controlF)]
print (paste("Done with ",x));
return (y)
}

# chromosome specific chip data filtering using ShortRead filter
extChr_Chip=function(x){
chrs=paste("chr",x,"$",sep="")
y=chipF[(chromosomeFilter(chrs))(chipF)]
print (paste("Done with ",x));
return (y)
}

# List of extracted chromosomes - Control
for (i in 1:19){
controlChr[[i]]=lapply(i,extChr_Control)
}

# List of extracted chromosomes - Chip
for (i in 1:19){
chipChr[[i]]=lapply(i,extChr_Chip)
}

# filtering X,Y and M chromosome
# control
controlChr[[20]]=controlF[which(controlF@chromosome=="chrX")]
controlChr[[21]]=controlF[which(controlF@chromosome=="chrY")]
controlChr[[22]]=controlF[which(controlF@chromosome=="chrM")]

# chip
chipChr[[20]]=chipF[which(chipF@chromosome=="chrX")]
chipChr[[21]]=chipF[which(chipF@chromosome=="chrY")]
chipChr[[22]]=chipF[which(chipF@chromosome=="chrM")]

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
chrs<-mapply(extractChr,control1level=controllevel,chip1level=chiplevel,MoreArgs=list(control1=control1,chip1=chip1))

# unlisting chip and control chromosome data to remove head list
controlChr=unlist(controlChr)
chipChr=unlist(chipChr)

# variable declarations
chip.norm.aln<-list(); control.norm.aln<-list()
chip.fragment.size<-c(); control.fragment.size<-c()
chip.cov<-c(); control.cov<-c()

# NORMALIZING LIBRARY SIZE FROM CHIP AND CONTROL
for (i in 1:length(chipChr)){
sizeControl<-length(controlChr[[i]]@id)
sizeChip<-length(chipChr[[i]]@id)
smallest.size<-min(sizeControl,sizeChip)
chip.norm.aln[[i]]<-sample(chipChr[[i]],smallest.size,replace=FALSE)
control.norm.aln[[i]]<-sample(controlChr[[i]],smallest.size,replace=FALSE)
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
for (i in 1:length(chipChr)){
chip.cov<-c(chip.cov,coverage(chip.norm.aln[[i]],width=seqlengths(Mmusculus),extend=as.integer(chip.extend.length[i])))
control.cov<-c(control.cov,coverage(control.norm.aln[[i]],width=seqlengths(Mmusculus),extend=as.integer(control.extend.length[i])))
}


# plotting function for coverage visualization of each chromosome

##########################################
# End of Code


