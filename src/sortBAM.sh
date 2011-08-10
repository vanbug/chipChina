#!bin/sh
echo "This script will sort all the bam files in the user defined folder using samtools"
for i in `ls *.bam`
do
echo "Sorting $i"
        samtools sort /projects/globalscratch/sukhi/Beijing/data/mapping/runII/bam/$i /projects/globalscratch/sukhi/Beijing/data/mapping/runII/sorted/$i.sort
done

