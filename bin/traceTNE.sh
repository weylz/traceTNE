#!/bin/bash
# File    : traceTNE.sh
# Time    : 2024/10/01 08:08:00
# Author  : Wenyong Zhu
# Version : 1.2.0
# Desc    : To trace the transcribed noncoding elements in sample(s) according to the defined criteria.


# Environment:
#     - Python 3.11.10;
#     - R 4.2.3 (library fitdistrplus 1.1-11);
#     - BEDTools suite v2.31.0;
#     - Jim Kent's executable programms: http://hgdownload.cse.ucsc.edu/admin/exe/.

# $@-List of all parameters
Options=$@
# $#-The number of arguments added to the Shell
Optnum=$#

# $0-The file name of the Shell itself ($0 $Options)
function usage(){
cat << EOF

Usage:    ${0##*\/} -I <...> -B <...> -S <...> [-L <...>] [-M <...>] [-V] [-H]

Options:
    -H    Print this help menu.
    -I    Label for input sample(s) (the name of the folder where the bigwig files are placed).
    -B    Bed file for blocklist (without a header requires three columns for chromosome, start, end).
    -S    Genome chromosome sizes (hg38.chrom.sizes or hg19.chrom.sizes).
    -L    Minimal length of transcribed noncoding elements (default: 150 bp).
    -M    Minimum value (included) of the transcriptional signal in the sample(s) (default: average transcriptional signal).
    -V    The version of this tool.
EOF
}

function scriptInfo(){
    if [[ $1 == "ERROR" ]]; then
        echo "-------------------------------------------------"
        echo ">>>>>>>>>>>>>>>>>>>>> ERROR <<<<<<<<<<<<<<<<<<<<<"
        echo "-------------------------------------------------"
    elif [[ $1 == "START" ]]; then
        echo "-------------------------------------------------"
        echo ">>>>>>>>>>>>>>>>>>>>> Start <<<<<<<<<<<<<<<<<<<<<"
        echo "-------------------------------------------------"
    elif [[ $1 == "END" ]]; then
        echo "-------------------------------------------------"
        echo ">>>>>>>>>>>>>>>>>>>>>> END <<<<<<<<<<<<<<<<<<<<<<"
        echo "-------------------------------------------------"
    elif [[ $1 == "VERSION" ]]; then
        echo -e "traceTNE v1.2.0\n"
    fi
}

# $?-Closing code of the last command to run (return value)
if [ $# = 0 ]; then usage && exit $?; fi

# ;;-Use the options in case and play the role of Terminator
while getopts "I:B:S:L:M:HV" OPTION
do
    case "$OPTION" in
        H ) usage && exit $?               ;;
        I ) sample_label="$OPTARG"         ;;
        B ) block_list="$OPTARG"           ;;
        S ) chrom_size="$OPTARG"           ;;
        L ) min_len="$OPTARG"              ;;
        M ) min_sig="$OPTARG"              ;;
        V ) scriptInfo VERSION && exit $?  ;;
        * ) scriptInfo ERROR && exit $?    ;;
    esac
done

[ -z "$sample_label" ] || [ ! -d "$sample_label"  ] && echo "Error: Label not found." && scriptInfo ERROR && exit $?
[ -z "$block_list" ] || [ ! -e "$block_list" ] && echo "Error: Blocklist not found." && scriptInfo ERROR && exit $?
[ -z "$chrom_size" ] || [ ! -e "$chrom_size" ] && echo "Error: File chromosome size not found." && scriptInfo ERROR && exit $?
[ -z "$min_len" ] && min_len=150
if [[ `echo "$min_len<=0" | bc` -eq 1 ]]; then echo "Error: Unacceptable element length." && scriptInfo ERROR && exit $?; fi
[ -z "$min_sig" ] || if [[ `echo "$min_sig<0" | bc` -eq 1 ]]; then echo "Error: Unacceptable transcriptional signaling threshold." && scriptInfo ERROR && exit $?; fi
#Exit if file already exists
[ -e "TNE."$sample_label".bed" ] && echo "File TNE."$sample_label".bed already exists." && exit $?

scriptInfo START && echo -e "\nInput: ${0##*\/} $@ \n"

echo "["`date`"] Scanning sample files in the "$sample_label" folder ..."
cd "$sample_label" && list_sample=($(ls -1 *.bw | sed 's/\.bw$//')) && num_sample=${#list_sample[@]}
[ $num_sample == 0 ] && echo -e "File(s) not found in the "$sample_label" folder." && scriptInfo ERROR && exit $?
echo -e "\tNumber of Sample(s) = "$num_sample
# if not exist, then move sample file & generate list_bigwig_$sample_label.txt
[ -e list_bigwig_"$sample_label".txt ] || for tmp in $(ls ./); do mkdir ${tmp%%.*} && mv $tmp ${tmp%%.*}; echo -e "$PWD/${tmp%%.*}/$tmp" >> list_bigwig_"$sample_label".txt; done

echo "["`date`"] Measuring average transcriptional signal in genomic regions ..."
CHROM_SIZE=`cat $chrom_size | awk '$1 ~ /^(chr[0-9]+|chrX|chrY)$/ {print}' | awk '{s+=$2}END{printf "%.0f\n", s}'`
if [ "$num_sample" == "1" ]
then
    bigWigToBedGraph $(cat list_bigwig_"$sample_label".txt) stdout | LC_ALL=C sort -S 20G -T /tmp -k1,1 -k2,2n | awk '$1 ~ /^(chr[0-9]+|chrX|chrY)$/ {print}' > trimmed_sig_"$sample_label".bg
else
    bigWigMerge -inList list_bigwig_"$sample_label".txt stdout | LC_ALL=C sort -S 20G -T /tmp -k1,1 -k2,2n | awk -v NUM=$num_sample '{FS=OFS="\t"; print $1,$2,$3,$4/NUM}' | awk '$1 ~ /^(chr[0-9]+|chrX|chrY)$/ {print}' > trimmed_sig_"$sample_label".bg
fi
avg_sig=`cat trimmed_sig_"$sample_label".bg | awk -v CHROM_SIZE=$CHROM_SIZE 'BEGIN{FS="\t"; S=0} {S+=$4*($3-$2)} END{print S/CHROM_SIZE}'`

echo "["`date`"] Detecting regions with higher RNAseq density than the average transcriptional siganl ..."
echo -e "\tAverage Transcriptional Signal = "$avg_sig
awk -v AVG_SIG=$avg_sig '{OFS="\t"; if($4>=AVG_SIG) print $1,$2,$3,$4}' trimmed_sig_"$sample_label".bg | mergeBed -c 4 -o max > TNE."$sample_label".1.tmp

echo "["`date`"] Extracting regions with transcriptional signal at significiant level (p ≤ 0.05) ..."
[ -e transcriptional_noise.txt ] || bedtools random -seed 1 -g $chrom_size -l 1 -n 1000000 | sortBed | intersectBed -a - -b $block_list -v -sorted | intersectBed -a trimmed_sig_"$sample_label".bg -b - -sorted -u | cut -f4 > transcriptional_noise.txt
fitRTS.R transcriptional_noise.txt
distp=`tail -n1 transcriptional_noise_pvalues.txt`
echo -e "\tSignificiant Level of Transcriptional Siganl = $distp"
awk -v DISTP=$distp '{OFS="\t"; if($4>=DISTP) print $1,$2,$3,$4}' TNE."$sample_label".1.tmp | mergeBed -d $min_len -c 4 -o max > TNE."$sample_label".2.tmp

echo "["`date`"] Removing regions that overlap the blocklist ..."
bedtools subtract -a TNE."$sample_label".2.tmp -b $block_list > TNE."$sample_label".3.tmp

echo "["`date`"] Restricting to a minimum length of "$min_len" bp ..."
awk -v MIN_LEN=$min_len '{OFS="\t"; if(($3-$2)>=MIN_LEN) print $1,$2,$3,$1"_"$2"_"$3}' TNE."$sample_label".3.tmp > TNE."$sample_label".4.tmp

echo "["`date`"] Creating random transcriptional regions and calculating their signals ..."
while read sample_bigwig
do
    count=$(grep -n "$sample_bigwig" list_bigwig_"$sample_label".txt | cut -d ':' -f 1); total=$(wc -l < list_bigwig_"$sample_label".txt); base="${sample_bigwig##*/}"; echo "["`date`"] ["$count"/"$total"]" "${base%.*}"
    bedtools shuffle -seed 1 -excl $block_list -i TNE."$sample_label".4.tmp -noOverlapping -g $chrom_size | awk -v OFS="\t" '$4=$1"_"$2"_"$3' | bigWigAverageOverBed $sample_bigwig stdin stdout | cut -f1,5 > $sample_bigwig.rdsig
    bigWigAverageOverBed $sample_bigwig TNE."$sample_label".4.tmp stdout | cut -f1,5 | sort -k1,1 > $sample_bigwig.avgsig
done < list_bigwig_"$sample_label".txt

echo "["`date`"] Calculating p value for TNE with binomial test (p ≤ 0.05) and FDR adjustment ..."
calSIG.R $sample_label list_bigwig_"$sample_label".txt

echo "["`date`"] Selecting TNE with FDR-adjusted p value ≤ 0.05 ..."
awk '{OFS="\t"; split($1,a,"_"); if($1~/^chr/) {if($5<=0.05) print a[1],a[2],a[3],$1}}' "TNE."$sample_label".adj_pvalues.txt" | sortBed -i stdin > TNE."$sample_label".bed

echo "["`date`"] ALL FINISHED."
scriptInfo END
