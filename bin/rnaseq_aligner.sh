#!/bin/bash
# File    : script_RNAseq_aligner.sh
# Time    : 2024/10/01 08:08:08
# Author  : Wenyong Zhu
# Version : 1.8.0
# Desc    : Note the valid directory locations and the value of 'sjdbOverhang' in the STAR index.


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
echo [ `date` ] "SCRIPT STARTING ..."
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
genome_fa="/path/to/xxx.fa"
genome_gtf="/path/to/xxx.gtf"
software_trimmomatic="/path/to/xxx.jar"
target_directory="$(pwd)/" # workspace for current project which contain raw data folder (00_rawData)
temp_directory="./STARtmp/" # temp directory for STAR (outTmpDir)
num_processor=$(nproc)
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Exit if the raw data folder does not exist.
[ ! -d 00_rawData ] && echo "Unable to locate the raw data folder. Please move all SRA files to folder (00_rawData) in the same directory as this script!" && exit $?
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# This script initiates the RNA-seq alignment from the SRA file.
echo [ `date` ] "Collecting sample information ..."
list=()
for filename in $(ls $target_directory/00_rawData)
do
	list+=($filename)
done
num_sample=${#list[*]}
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
echo [ `date` ] "fastq-dump / fasterq-dump starting ..."
cd $target_directory
[ ! -e 01_fastq ] && mkdir 01_fastq
cd 01_fastq
tmp_fifofile="/tmp/$$.fifo"; mkfifo $tmp_fifofile; exec 6<>$tmp_fifofile; rm $tmp_fifofile
thread_num=5; for ((i=0;i<${thread_num};i++));do echo; done >&6
num=0
for filename in ${list[*]}
do
	num=`expr $num + 1`
	tmp_list=`ls $target_directory/01_fastq`
	if [[ "${tmp_list[@]}"  =~ "$filename" ]]
	then
		echo [ `date` ] $num/$num_sample $filename "fastq-dump / fasterq-dump | gzip / pigz exists --> Skip."
	else
		read -u6
		{
			echo [ `date` ] $num/$num_sample $filename "fastq-dump / fasterq-dump | gzip / pigz ..."
			fasterq-dump --threads $((num_processor / 2)) --mem 30G --bufsize 20G --curcache 20G --split-3 --outdir $target_directory/01_fastq/ $target_directory/00_rawData/$filename
			for i in `seq 1 2`
			do
			pigz -p $((num_processor / 2)) -q $target_directory/01_fastq/$filename"_"$i".fastq" &
			done
			echo [ `date` ] $num/$num_sample $filename "fastq-dump / fasterq-dump | gzip / pigz finished successfully."
			wait
			echo >&6
		} &
	fi
done
wait
exec 6>&-
echo [ `date` ] "fastq-dump / fasterq-dump done."
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
echo [ `date` ] "Checking the format of the raw data ..."
cd $target_directory/01_fastq
num=0
for filename in ${list[*]}
do
    num=`expr $num + 1`
	if [[ $filename == *.fastq.gz ]]
	then
		echo [ `date` ] $num/$num_sample $filename " PASS."
	fi
done
len_sequence=`zcat "$target_directory"/01_fastq/"$list"_1.fastq.gz | awk 'NR%4==2 {print length($1)}'| head -n1`
echo [ `date` ] "All raw data formatting checked."
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
echo [ `date` ] "fastqc starting ..."
cd $target_directory
[ ! -e 02_fastqc ] && mkdir 02_fastqc
cd 02_fastqc
tmp_fifofile="/tmp/$$.fifo"; mkfifo $tmp_fifofile; exec 6<>$tmp_fifofile; rm $tmp_fifofile
thread_num=5; for ((i=0;i<${thread_num};i++));do echo; done >&6
num=0
for filename in ${list[*]}
do
	num=`expr $num + 1`
	tmp_list=`ls $target_directory/02_fastqc`
	if [[ "${tmp_list[@]}"  =~ "$filename" ]]
	then
		echo [ `date` ] $num/$num_sample $filename "fastqc exists --> Skip."
	else
		read -u6
		{
			echo [ `date` ] $num/$num_sample $filename "fastqc ..."
			mkdir $filename
			for i in `seq 1 2`
			do
			fastqc -t $((num_processor / 2)) -q -o "./"$filename"/" --nogroup $target_directory"/01_fastq/"$filename"_"$i".fastq.gz" &
			done
			echo [ `date` ] $num/$num_sample $filename "fastqc finished successfully."
			wait
			echo >&6
		} &
	fi
done
wait
exec 6>&-
echo [ `date` ] "fastqc done."
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
echo [ `date` ] "trimmomatic starting ..."
cd $target_directory
[ ! -e 03_trimmomatic ] && mkdir 03_trimmomatic
cd 03_trimmomatic
tmp_fifofile="/tmp/$$.fifo"; mkfifo $tmp_fifofile; exec 6<>$tmp_fifofile; rm $tmp_fifofile
thread_num=5; for ((i=0;i<${thread_num};i++));do echo; done >&6
num=0
for filename in ${list[*]}
do
	num=`expr $num + 1`
	tmp_list=`ls $target_directory/03_trimmomatic`
	if [[ "${tmp_list[@]}"  =~ "$filename" ]]
	then
		echo [ `date` ] $num/$num_sample $filename "trimmomatic exists --> Skip."
	else
		read -u6
		{
			echo [ `date` ] $num/$num_sample $filename "trimmomatic ..."
			java -jar $software_trimmomatic PE -threads $((num_processor / 2)) -phred33 -quiet $target_directory"/01_fastq/"$filename"_1.fastq.gz" $target_directory"/01_fastq/"$filename"_2.fastq.gz" $target_directory"/03_trimmomatic/"$filename"_1_paired.fastq.gz" $target_directory"/03_trimmomatic/"$filename"_1_unpaired.fastq.gz" $target_directory"/03_trimmomatic/"$filename"_2_paired.fastq.gz" $target_directory"/03_trimmomatic/"$filename"_2_unpaired.fastq.gz" ILLUMINACLIP:$software_trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
			echo [ `date` ] $num/$num_sample $filename "trimmomatic finished successfully."
			echo >&6
		} &
	fi
done
wait
exec 6>&-
echo [ `date` ] "trimmomatic done."
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
echo [ `date` ] "STAR starting ..."
cd $target_directory
[ ! -e 04_STAR ] && mkdir 04_STAR
cd 04_STAR
if [ -e star_index ]
then
	echo [ `date` ] "STAR index exists --> Skip."
else
	echo [ `date` ] "STAR index generating ..."
	mkdir star_index && cd star_index
	STAR --runThreadN $num_processor --runMode genomeGenerate --genomeDir $target_directory/04_STAR/star_index/ --genomeFastaFiles $genome_fa --sjdbGTFfile $genome_gtf --sjdbOverhang $((len_sequence - 1))
fi
cd $target_directory/04_STAR
if [ ! -e $temp_directory ]; then mkdir $temp_directory; else rm -rf $temp_directory/*; fi
ulimit -n 65535
tmp_fifofile="/tmp/$$.fifo"; mkfifo $tmp_fifofile; exec 6<>$tmp_fifofile; rm $tmp_fifofile
thread_num=5; for ((i=0;i<${thread_num};i++));do echo; done >&6
num=0
for filename in ${list[*]}
do
	num=`expr $num + 1`
	tmp_list=`ls $target_directory/04_STAR`
	if [[ "${tmp_list[@]}"  =~ "$filename" ]]
	then
		echo [ `date` ] $num/$num_sample $filename "STAR exists --> Skip."
	else
		read -u6
		{
			echo [ `date` ] $num/$num_sample $filename "STAR ..."
			mkdir $filename
			STAR --runThreadN $((num_processor / 2)) --genomeDir $target_directory/04_STAR/star_index --readFilesIn $target_directory"/03_trimmomatic/"$filename"_1_paired.fastq.gz" $target_directory"/03_trimmomatic/"$filename"_2_paired.fastq.gz" --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $target_directory"/04_STAR/"$filename"/"$filename"_" --outTmpDir $temp_directory/$filename --readFilesCommand zcat
			echo [ `date` ] $num/$num_sample $filename "STAR finished successfully."
			echo >&6
		} &
	fi
done
wait
exec 6>&-
rm -rf $temp_directory
echo [ `date` ] "STAR done."
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
echo [ `date` ] "samtools starting ..."
cd $target_directory
[ ! -e 05_samtools_sort ] && mkdir 05_samtools_sort
cd 05_samtools_sort
tmp_fifofile="/tmp/$$.fifo"; mkfifo $tmp_fifofile; exec 6<>$tmp_fifofile; rm $tmp_fifofile
thread_num=5; for ((i=0;i<${thread_num};i++));do echo; done >&6
num=0
for filename in ${list[*]}
do
	num=`expr $num + 1`
	tmp_list=`ls $target_directory/05_samtools_sort`
	if [[ "${tmp_list[@]}"  =~ "$filename" ]]
	then
		echo [ `date` ] $num/$num_sample $filename "samtools sort & index exists --> Skip."
	else
		read -u6
		{
			echo [ `date` ] $num/$num_sample $filename "samtools sort & index ..."
			cd $target_directory/05_samtools_sort && mkdir $filename && cd $filename
			samtools sort -@ $((num_processor / 2)) $target_directory"/04_STAR/"$filename"/"$filename"_Aligned.sortedByCoord.out.bam" -o $target_directory"/05_samtools_sort/"$filename"/"$filename"_Aligned.sortedByCoord.out.sort.bam"
			samtools index -@ $((num_processor / 2)) -b $target_directory"/05_samtools_sort/"$filename"/"$filename"_Aligned.sortedByCoord.out.sort.bam"
			echo [ `date` ] $num/$num_sample $filename "samtools sort & index finished successfully."
			echo >&6
		} &
	fi
done
wait
exec 6>&-
echo [ `date` ] "samtools done."
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
echo [ `date` ] "bamCoverage starting ..."
cd $target_directory
[ ! -e 06_bamCoverage ] && mkdir 06_bamCoverage
cd 06_bamCoverage
tmp_fifofile="/tmp/$$.fifo"; mkfifo $tmp_fifofile; exec 6<>$tmp_fifofile; rm $tmp_fifofile
thread_num=5; for ((i=0;i<${thread_num};i++));do echo; done >&6
num=0
for filename in ${list[*]}
do
	num=`expr $num + 1`
	tmp_list=`ls $target_directory/06_bamCoverage`
	if [[ "${tmp_list[@]}"  =~ "$filename" ]]
	then
		echo [ `date` ] $num/$num_sample $filename "bamCoverage exists --> Skip."
	else
		read -u6
		{
			echo [ `date` ] $num/$num_sample $filename "bamCoverage ..."
			bamCoverage -p $((num_processor / 2)) -b $target_directory"/05_samtools_sort/"$filename"/"$filename"_Aligned.sortedByCoord.out.sort.bam" -o $target_directory"/06_bamCoverage/"$filename".bw" --normalizeUsing BPM --binSize 10
			echo [ `date` ] $num/$num_sample $filename "bamCoverage finished successfully."
			echo >&6
		} &
	fi
done
wait
exec 6>&-
echo [ `date` ] "bamCoverage done."
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
echo [ `date` ] "ALL FINISHED SUCCESSFULLY."
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
