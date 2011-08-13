# shell script for the sam single end file production
#!bin/sh
files=(`ls /projects/globalscratch/sukhi/beijing/data/fastq/*.clean.fq`)
sai=(`ls *.sai`)
size=${#files[@]}

# starting samse (sam production for single end reads)
for ((i=0; i<size; i++));
do
  	echo "Producing ${files[i]}.sam from ${sai[i]}"
	nice bwa samse /biodata/biodb/ABG/genomes/mouse/mm9/mm9 ${sai[i]} ${files[i]} > /projects/globalscratch/sukhi/beijing/data/mapping/bwa/sam/$i.sam
	#echo ${files[i]}
	#echo ${sai[i]}
done
