#!bin/sh 
echo "This script will convert all the sam files to bam files in the user defined folder in this script 
using samtools"
for i in `ls *.sam`
do
echo "Converting $i to bam (machine readable & compressed)"
	samtools view -bS -o /projects/globalscratch/sukhi/Beijing/data/mapping/runII/bam/$i.bam /projects/globalscratch/sukhi/Beijing/data/mapping/runII/sam/$i
done
