###file_path
input_path=${all_path}/Raw_data/${datainfo}
file_path=`pwd`
raw_dir=${file_path}/00_rawdata
logs_dir=${file_path}/logs
fastqc_dir=${file_path}/01_fastqc
trimmedFastq_dir=${file_path}/00_trimmeiddata
trimmedFastq_log_dir=${file_path}/logs/trimmeddata
trimmed_fastqc_dir=${file_path}/02_trimmed_fastqc
align_exp_dir=${file_path}/03_bam_hg38
alignexp_log_dir=${file_path}/logs/align_hg38
exp_bam_rmdup=${file_path}/03_bam_hg38_rmdup
rmdup_exp_log=${file_path}/logs/rmdup_state_hg38
align_spike_dir=${file_path}/03_spikebam_mm10
alignspike_log_dir=${file_path}/logs/align_mm10
spike_bam_rmdup=${file_path}/03_spikebam_mm10_rmdup
rmdup_spike_log=${file_path}/logs/rmdup_state_mm10
bw_dir=${file_path}/04_bw_rmdup
sampleinfo=${file_path}/sampleinfo.txt
trimmed_5prime=${file_path}/02_1_trimmed_5prime
trimmed_5prime_log=${file_path}/logs/02_1_trimmed_5prime_log
trimmed_5prime_fastq=${file_path}/02_1_trimmed_5prime_fastq

#reference genome
GENOME_EXP="${file_path}/hg38"
GENOME_SPIKE="${file_path}/mm10"
SPIKE_PREFIX="Mus"

MAPQ=10

#step 0
####check files 

#step 1.1
####change file name#####
mkdir -p ${raw_dir}
ln -s `ls ${input_path}/*/*.gz` ${raw_dir}/ >/dev/null 2>&1
#
cat $sampleinfo| while read id;
do
arr=($id)
sample1=${arr[0]}
sample2=${arr[1]}
echo "sample1:$sample1 sample2:$sample2"
fq1=$(ls ${raw_dir}/*raw_1.fq.gz |grep "$sample1") 
fq2=$(ls ${raw_dir}/*raw_2.fq.gz |grep "$sample1") 
mv $fq1 ${raw_dir}/${sample2}_R1.fastq.gz
mv $fq2 ${raw_dir}/${sample2}_R2.fastq.gz
done
#
#
##step 1.2
#####fastqc of raw data ####
mkdir -p ${fastqc_dir}
#
in_path=${raw_dir}
out_path=${fastqc_dir}
#
for FILE in `ls ${in_path}/*.gz`
do 
	if [ ! -s ${out_path}/"$(basename ${FILE/.fastq.gz/_fastqc.zip})" ]
	then
	   fastqc $FILE -t 2 -o ${out_path}/ &
	fi
done
wait
#
##step 1.3
##check fq
cd ${fastqc_dir}
ls *fastqc.zip|sed 's/_R[1-2]_fastqc.zip//g'|sort|uniq -c
#
#### merge reports of fastqc
     multiqc ${fastqc_dir}/ -n rawdata_multiqc -o ${fastqc_dir}/
#
#
#step 2.1
### Trimming adapters  (trim_galore)
####remember to change --length for proseq !!!!!! #####
#
mkdir -p ${trimmedFastq_dir}
mkdir -p ${trimmedFastq_log_dir}
#
for fq1 in `ls ${raw_dir}/*R1.fastq.gz`

do
fq2=${fq1/R1.fastq.gz/R2.fastq.gz}
    if [ ! -s ${trimmedFastq_dir}/"$(basename ${fq1/.fastq.gz/_val_1.fq.gz})" ]
    then
    trim_galore --paired -o ${trimmedFastq_dir} $fq1 $fq2  \
    > ${trimmedFastq_log_dir}/"$(basename ${fq1/_R1.fq.gz/_trimmed.log})" 2>&1 & 
    fi
done
wait
#
##step 2.4
######qc for trimmed data####
mkdir -p ${trimmed_fastqc_dir}
#
for FILE in `ls ${trimmedFastq_dir}/*fq.gz`
do
if [ ! -s ${trimmed_fastqc_dir}/"$(basename ${FILE/.fq.gz/_fastqc.html})" ]
then
 fastqc $FILE -t 2 -o ${trimmed_fastqc_dir} &
fi
done
wait
#### merge reports of fastqc
multiqc ${trimmed_fastqc_dir}/ -n trimmed_multiqc -o ${trimmed_fastqc_dir}/






#step 3.1
###Aligning to experimental genome remove chrM reads#####
mkdir ${align_exp_dir}
mkdir ${alignexp_log_dir}

for PAIR in $(ls ${trimmedFastq_dir} | sed 's/_R[1-2].*//' | sort | uniq )
do
if [ ! -s "${align_exp_dir}/${PAIR}_hg38.bam" ]
then
    echo "aligning ${PAIR} to experimental genome"
    (bowtie2 \
    -p 25 -t -q -N 1 -L 25 -X 2000 --no-mixed --no-discordant \
    -x "$GENOME_EXP" \
    -1 "${trimmedFastq_dir}/${PAIR}_R1_val_1.fq.gz" \
    -2 "${trimmedFastq_dir}/${PAIR}_R2_val_2.fq.gz" \
    2> ${alignexp_log_dir}/${PAIR}_align.log) |
    samtools view -bS -f 2 -q ${MAPQ} |
    samtools sort -@ 20 -o ${align_exp_dir}/${PAIR}_hg38.bam
    samtools index ${align_exp_dir}/${PAIR}_hg38.bam
fi
done

#step 3.2
### Aligning to spike in genome to get normalization factors ###
mkdir -p ${align_spike_dir}
mkdir -p ${alignspike_log_dir}

for PAIR in $(ls ${trimmedFastq_dir} | sed 's/_R[1-2].*//' | sort | uniq )
do
if [ ! -s "${align_spike_dir}/${PAIR}_mm10.bam" ]
    then
    echo "aligning ${PAIR} to spike-in genome"
    (bowtie2 \
    -x "$GENOME_SPIKE" \
    -p 25 -t -q -N 1 -L 25 -X 2000 --no-mixed --no-discordant \
    -1 "${trimmedFastq_dir}/${PAIR}_R1_val_1.fq.gz" \
    -2 "${trimmedFastq_dir}/${PAIR}_R2_val_2.fq.gz" \
    2> ${alignspike_log_dir}/${PAIR}_mm10Align.log) |
    samtools view -bS -f 2 -q ${MAPQ} |
    samtools sort -@ 20 -o ${align_spike_dir}/${PAIR}_mm10.bam
    samtools index ${align_spike_dir}/${PAIR}_mm10.bam
fi
done

### step 3.3
## remove dulipate/mito and shift
###remove duplicates of exp genome
mkdir -p ${exp_bam_rmdup}
mkdir -p ${rmdup_exp_log}

cd ${exp_bam_rmdup}
ls ${align_exp_dir}/*.bam|while read id;
do
sample=$(basename ${id/.bam/})
if [ ! -s "${exp_bam_rmdup}/${sample}.rmdup.bam" ]
    then
       picard  MarkDuplicates -REMOVE_DUPLICATES True \
        -I $id \
        -O ${exp_bam_rmdup}/${sample}.rmdup.bam \
        -M ${exp_bam_rmdup}/${sample}.rmdup.metrics
        samtools index ${exp_bam_rmdup}/${sample}.rmdup.bam
        samtools flagstat ${exp_bam_rmdup}/${sample}.rmdup.bam > ${rmdup_exp_log}/${sample}.rmdup.stat
    fi
	samtools view -h ${exp_bam_rmdup}/${sample}.rmdup.bam | grep -v chrM | samtools sort -O bam  -@ 25 -o - > ${exp_bam_rmdup}/${sample}.rmdup_rmChrM.bam
	samtools index ${exp_bam_rmdup}/${sample}.rmdup_rmChrM.bam
	alignmentSieve --numberOfProcessors 25 --ATACshift -b ${exp_bam_rmdup}/${sample}.rmdup_rmChrM.bam -o ${exp_bam_rmdup}/${sample}.rmdup_rmChrM_Shift.bam 
	samtools sort ${exp_bam_rmdup}/${sample}.rmdup_rmChrM_Shift.bam -O bam -@ 25 -o ${exp_bam_rmdup}/${sample}.rmdup_rmChrM_Shift.sorted.bam
	samtools index ${exp_bam_rmdup}/${sample}.rmdup_rmChrM_Shift.sorted.bam
	samtools flagstat ${exp_bam_rmdup}/${sample}.rmdup_rmChrM_Shift.sorted.bam > ${rmdup_exp_log}/${sample}.rmdup_rmChrM_Shift.sorted.stat
done
#
# step 3.4
#remove duplicates of spike-in genome
mkdir -p ${spike_bam_rmdup}
mkdir -p ${rmdup_spike_log}

cd ${spike_bam_rmdup}
ls ${align_spike_dir}/*.bam|while read id;
do
sample=$(basename ${id/.bam/})
if [ ! -s "${spike_bam_rmdup}/${sample}.rmdup.bam" ]
    then
        picard  MarkDuplicates -REMOVE_DUPLICATES True \
        -I $id \
        -O ${spike_bam_rmdup}/${sample}.rmdup.bam \
        -M ${spike_bam_rmdup}/${sample}.rmdup.metrics
        samtools index ${spike_bam_rmdup}/${sample}.rmdup.bam
        samtools flagstat ${spike_bam_rmdup}/${sample}.rmdup.bam > ${rmdup_spike_log}/${sample}.rmdup.stat
    fi
done

#step 3.5
### calculate normalization factors ###
hg38_path=${rmdup_exp_log}
mm10_path=${rmdup_spike_log}
align_path=${alignexp_log_dir}

echo -e "sample\tALLREADS\thg38_READS\thg38_mapping_RATIO\thg38_qc_READS\thg38_qc_RATIO\thg38_rmChrM_READS\thg38_rmChrM_RATIO\tmtDNA_READS\tmtDNA_RATIO\tmm10_qc_READS\tmm10_qc_RATIO_intotal\tmm10_qc_RATIO_inqc\tSCALEFACTOR\t" >${logs_dir}/scalefactor.txt
ls $align_path|while read file;
do
    sample=${file/_align.log/}
    ALLREADS=$(cat ${align_path}/$file|grep "were paired; of these:$"|cut -d "(" -f 1|awk '{print $1*2}')
    hg38_READS=$(cat ${align_path}/$sample"_align.log"| sed 's/%//g' | awk '{printf $0"\t"}'  |cut -f 8,9 | sed 's/\t/\n/g' | awk '{print $1}' | awk '{printf $0"\t"}'|awk '{print 2*($1+$2)}')
    hg38_mapping_RATIO=$(cat ${align_path}/$file|grep "overall alignment rate"|cut -d "%" -f 1)
    hg38_qc_READS=$(cat $hg38_path/${sample}_hg38.rmdup.stat|grep "total (QC-passed reads"|cut -d " " -f 1)
    hg38_qc_RATIO=`printf "%.2f\n" $(echo "${hg38_qc_READS}/${ALLREADS}*100"|bc -l)`
    hg38_rmChrM_READS=$(cat $hg38_path/${sample}"_hg38.rmdup_rmChrM_Shift.sorted.stat"|grep "total (QC-passed reads"|cut -d " " -f 1)
    hg38_rmChrM_RATIO=`printf "%.2f\n" $(echo "${hg38_rmChrM_READS}/${ALLREADS}*100"|bc -l)`
    mtDNA_READS=$(samtools idxstats ${align_exp_dir}/${sample}"_hg38.bam" | grep 'chrM' | awk '{SUM += $3} END {print SUM}')
    mtDNA_RATIO=`printf "%.2f\n" $(echo "${mtDNA_READS}/${ALLREADS}*100"| bc -l)`
    mm10_qc_READS=$(cat ${mm10_path}/${sample}_mm10.rmdup.stat|grep "total (QC-passed reads"|cut -d "+" -f 1)
    QC_reads=$(echo "${mm10_qc_READS}+${hg38_qc_READS}"|bc )
    mm10_qc_RATIO_intotal=`printf "%.2f\n" $(echo "${mm10_qc_READS}/${ALLREADS}*100"|bc -l)`
    mm10_qc_RATIO_inqc=`printf "%.2f\n" $(echo "${mm10_qc_READS}/${QC_reads}*100"|bc -l)`
    SCALEFACTOR=$(echo "1000000/${mm10_qc_READS}"|bc -l)
    echo -e $sample"\t"$ALLREADS"\t"$hg38_READS"\t"$hg38_mapping_RATIO"\t"$hg38_qc_READS"\t"$hg38_qc_RATIO"\t"$hg38_rmChrM_READS"\t"$hg38_rmChrM_RATIO"\t"$mtDNA_READS"\t"$mtDNA_RATIO"\t"$mm10_qc_READS"\t"$mm10_qc_RATIO_intotal"\t"$mm10_qc_RATIO_inqc"\t"$SCALEFACTOR >> ${logs_dir}/scalefactor.txt
done

echo "scale factor is done"

#step 4.1
### Making CPM-normalized bigWig files with reads adjusted by CPM ###
mkdir -p ${bw_dir}

cat  ${logs_dir}/scalefactor.txt | sed '1d' | while read id;
do
arr=($id)
sample=${arr[0]}

bam_file=${exp_bam_rmdup}/${sample}_hg38.rmdup_rmChrM_Shift.sorted.bam
    if [ ! -s "${bw_dir}/${sample}_CPM.bw" ]
    then
        bamCoverage \
        --bam ${bam_file} \
        --blackListFileName ${file_path}/hg38-blacklist.v2.bed \
        --outFileName ${bw_dir}/${sample}"_CPM.bw" \
        --binSize 1 \
        --scaleFactor 1 \
        --numberOfProcessors 23 \
        --normalizeUsing CPM
    fi
done

#step 4.2
### Making CPM-normalized bigWig files with reads adjusted by spike-in ###

cat  ${logs_dir}/scalefactor.txt | sed '1d' | while read id;
do
arr=($id)
sample=${arr[0]}
scalefactor=${arr[13]}
bam_file=${exp_bam_rmdup}/${sample}_hg38.rmdup_rmChrM_Shift.sorted.bam
    if [ ! -s "${bw_dir}/${sample}_spikein.bw" ]
    then
        bamCoverage \
        --bam ${bam_file} \
        --blackListFileName ${file_path}/hg38-blacklist.v2.bed \
        --outFileName ${bw_dir}/${sample}_spikein.bw \
        --binSize 1 \
        --scaleFactor $scalefactor \
        --numberOfProcessors 23 \
        --normalizeUsing None
    fi
done

#step 5.1
#### Call Peaks without input ###
peak_dir=${file_path}/05_peaks
peak_log=${file_path}/logs/peaks
mkdir -p ${peak_dir}/peak_alone
ls ${exp_bam_rmdup}/*_hg38.rmdup_rmChrM_Shift.sorted.bam | while read id;
do
    sample=`basename $id | sed 's/_hg38.rmdup_rmChrM_Shift.sorted.bam//g'`
    bam_file=${exp_bam_rmdup}/${sample}_hg38.rmdup_rmChrM_Shift.sorted.bam
    blacklist="${file_path}/hg38-blacklist.v2.bed"
    mkdir -p ${peak_log}/${sample}
    if [ ! -f ${peak_dir}/peak_alone/${sample}_peaks.narrowPeak ]
    then
        macs2 callpeak -t ${bam_file} -f BAMPE -n ${peak_dir}/peak_alone/$sample --nomodel --shift -75 --extsize 150  -g mm --keep-dup all --call-summits 2>${peak_log}/${sample}/${sample}_alone.log
        bedtools intersect -a ${peak_dir}/peak_alone/${sample}_peaks.narrowPeak -b $blacklist -f 0.25 -v > ${peak_dir}/peak_alone/${sample}_peaks.final.narrowPeak &
    fi
done
wait

