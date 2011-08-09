###################################################################################################################################
# R script for analysis of Chip Seq data from the Illumina for bowtie aligned data
# Author : Sukhdeep Singh
# Organization : Max Planck Dresden
###################################################################################################################################

library('ShortRead');
# declaring variables
chrs<-c();chipChr<-list();controlChr<-list();

# reading control and chip data
control<-readAligned('/projects/globalscratch/sukhi/Beijing/data/mapping/runII/sorted/GFP_Ab_Chip.clean.fq.sam.bam.sort.bam',type="BAM");
chip<-readAligned('/projects/globalscratch/sukhi/Beijing/data/mapping/runII/sorted/MEN1_Ab_Chip.clean.fq.sam.bam.sort.bam',type="BAM");

# printing info
print (control);
print (chip);

# uniquely aligned data, dropping NA's
controlP=control[which(is.na(control@strand)==F)]
chipP=chip[which(is.na(chip@strand)==F)]

# chromosome specific control data filtering using ShortRead filter - strange bug for chr-1 (picks all with 1 and doublets-10,11 etc)
extChr_Control=function(x){
chrs=paste("chr",x,"$",sep="")
y=controlP[(chromosomeFilter(chrs))(controlP)]
print (paste("Done with ",x));
return (y)
}

# chromosome specific chip data filtering using ShortRead filter
extChr_Chip=function(x){
chrs=paste("chr",x,"$",sep="")
y=chipP[(chromosomeFilter(chrs))(chipP)]
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
controlChr[[20]]=controlP[which(controlP@chromosome=="chrX")]
controlChr[[21]]=controlP[which(controlP@chromosome=="chrY")]
controlChr[[22]]=controlP[which(controlP@chromosome=="chrM")]

# chip
chipChr[[20]]=chipP[which(chipP@chromosome=="chrX")]
chipChr[[21]]=chipP[which(chipP@chromosome=="chrY")]
chipChr[[22]]=chipP[which(chipP@chromosome=="chrM")]


##########################################
# End of Code

