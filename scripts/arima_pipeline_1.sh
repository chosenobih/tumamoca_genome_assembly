#! /bin/bash
### Run Arima HiC mapping. Change parameters below and run with sh Arima_mapping_pipeline.sh

##############################################
# ARIMA GENOMICS MAPPING PIPELINE 07/26/2023 #
##############################################

# Below find the commands used to map HiC data.

# Replace the variables at the top with the correct paths for the locations of files/programs on your system.

# This bash script will map one paired end HiC dataset (read1 & read2 FASTQs). Feel to modify and multiplex as you see fit to work with your volume of samples and system.

##########################################
# Commands #
##########################################
SAMPLE_LIST=$1
if [ -z "${SAMPLE_LIST}" ]
then
    echo "This script needs the SAMPLE lists as the first parameter"
    echo "the SAMPLE info list is a CSV with 6 columns"
    echo "the 1st column is sample name which will be the prefix for file names"
    echo "the 2nd column is HiFi reads including path"
    echo "the 3rd column is Hi-C read pair 1 including path"
    echo "the 4th column is Hi-C read pair 2 including path"
    echo "Expected ploidy of sample"
    echo "Expected chromosome number"
    exit
else
    echo "Processing SAMPLEs from ${SAMPLE_LIST}"
fi
 
#---- PROCESSING
 
for SAMPLEINFO in $(cat $SAMPLE_LIST);
 
do
  NAME=$(echo ${SAMPLEINFO} | awk -F "\"*,\"*" '{print $1}')
  HiFi=$(echo ${SAMPLEINFO} | awk -F "\"*,\"*" '{print $2}')
  HC1=$(echo ${SAMPLEINFO} | awk -F "\"*,\"*" '{print $3}')
  HC2=$(echo ${SAMPLEINFO} | awk -F "\"*,\"*" '{print $4}')
  Ploidy=$(echo ${SAMPLEINFO} | awk -F "\"*,\"*" '{print $5}')
  Chr=$(echo ${SAMPLEINFO} | awk -F "\"*,\"*" '{print $6}')
  
ASSEMBLY="/data2/service_center_analyses/tumamoca_genome/raw_data/tumamoca_assembly_hifiasm.combined_racon.fa"
OUTPUT_NAME="${NAME}_Phase"

SRA=$NAME
LABEL=$NAME_Arima_Pipeline
BWA='bwa'
SAMTOOLS='samtools'
IN_DIR='Arima/trimmed'
REF=$ASSEMBLY
FAIDX=$REF.fai
PREFIX=$NAME
RAW_DIR='Arima'
FILT_DIR='Arima/filter'
FILTER='/data2/castri/Arima/filter_five_end.pl'
COMBINER='/data2/castri/Arima/two_read_bam_combiner.pl'
STATS='/data2/castri/Arima/get_stats.pl'
PICARD='/data2/castri/picard/build/libs/picard.jar'
TMP_DIR='Arima/temp'
TMP_DIR2='Arima/temp2'
TMP_DIR3='Arima/temp3'
PAIR_DIR='Arima/pair'
REP_DIR='Arima_hap12/rep'
REP_LABEL=${LABEL}_rep1
MERGE_DIR='Arima/merge'
MAPQ_FILTER=10
CPU=100
JAVA_TMP='tmp'

echo "### Step 0: Check output directories' existence & create them as needed"
[ -d $IN_DIR ] || mkdir -p $IN_DIR
[ -d $RAW_DIR ] || mkdir -p $RAW_DIR
[ -d $FILT_DIR ] || mkdir -p $FILT_DIR
[ -d $TMP_DIR ] || mkdir -p $TMP_DIR
[ -d $TMP_DIR2 ] || mkdir -p $TMP_DIR2
[ -d $TMP_DIR3 ] || mkdir -p $TMP_DIR3
[ -d $PAIR_DIR ] || mkdir -p $PAIR_DIR
[ -d $REP_DIR ] || mkdir -p $REP_DIR
[ -d $MERGE_DIR ] || mkdir -p $MERGE_DIR
[ -d $JAVA_TMP ] || mkdir -p $JAVA_TMP

echo "### Step 0: Index reference" # Run only once! Skip this step if you have already generated BWA index files
$BWA index -a bwtsw $REF

echo "Trim Hi-C reads" # Run only once
zcat $HC1| awk '{ if(NR%2==0) {print substr($1,6)} else {print} }' | gzip > $IN_DIR/${SRA}_1.fastq.gz
zcat $HC2| awk '{ if(NR%2==0) {print substr($1,6)} else {print} }' | gzip > $IN_DIR/${SRA}_2.fastq.gz

echo "### Step 1.A: FASTQ to BAM (1st)"
$BWA mem -t $CPU $REF $IN_DIR/${SRA}_1.fastq.gz > $RAW_DIR/${SRA}_1.sam
$SAMTOOLS view -@ $CPU -Sb $RAW_DIR/${SRA}_1.sam > $RAW_DIR/${SRA}_1.bam

echo "### Step 1.B: FASTQ to BAM (2nd)"
$BWA mem -t $CPU $REF $IN_DIR/${SRA}_2.fastq.gz  > $RAW_DIR/${SRA}_2.sam
$SAMTOOLS view -@ $CPU -Sb $RAW_DIR/${SRA}_2.sam > $RAW_DIR/${SRA}_2.bam

echo "### Step 2.A: Filter 5' end (1st)"
$SAMTOOLS view -h $RAW_DIR/${SRA}_1.bam | perl $FILTER | $SAMTOOLS view -Sb - > $FILT_DIR/${SRA}_1.bam

echo "### Step 2.B: Filter 5' end (2nd)"
$SAMTOOLS view -h $RAW_DIR/${SRA}_2.bam | perl $FILTER | $SAMTOOLS view -Sb - > $FILT_DIR/${SRA}_2.bam

echo "### Step 3A: Pair reads & mapping quality filter"
perl $COMBINER $FILT_DIR/${SRA}_1.bam $FILT_DIR/${SRA}_2.bam $SAMTOOLS $MAPQ_FILTER > $TMP_DIR/$SRA.bam

$SAMTOOLS view -bS -t $FAIDX $TMP_DIR/$SRA.bam > $TMP_DIR2/$SRA.bam

$SAMTOOLS sort -@ $CPU -o  $TMP_DIR3/$SRA.bam $TMP_DIR2/$SRA.bam

perl $COMBINER $FILT_DIR/${SRA}_1.bam $FILT_DIR/${SRA}_2.bam $SAMTOOLS $MAPQ_FILTER | $SAMTOOLS view -bS -t $FAIDX - | $SAMTOOLS sort -@ $CPU -o $TMP_DIR/$SRA.bam -

module load java

echo "### Step 3.B: Add read group"
#java -Xmx4G -Djava.io.tmpdir=temp/ -jar $PICARD AddOrReplaceReadGroups INPUT=$TMP_DIR3/${SRA}.bam OUTPUT=$PAIR_DIR/$SRA.bam ID=$SRA LB=$SRA SM=$LABEL PL=ILLUMINA PU=none

java -Xmx4G -Djava.io.tmpdir=tmp/ -jar $PICARD AddOrReplaceReadGroups INPUT=$TMP_DIR3/${SRA}.bam OUTPUT=$PAIR_DIR/$SRA.bam ID=$SRA LB=$SRA RGPL=ILLUMINA PU=none RGSM=20

###############################################################################################################################################################
###                                           How to Accommodate Technical Replicates                                                                       ###
### This pipeline is currently built for processing a single sample with one read1 and read2 FASTQ file.                                                    ###
### Technical replicates (eg. one library split across multiple lanes) should be merged before running the MarkDuplicates command.                          ###
### If this step is run, the names and locations of input files to subsequent steps will need to be modified in order for subsequent steps to run correctly.###
### The code below is an example of how to merge technical replicates.                                                                                      ###
###############################################################################################################################################################
#	REP_NUM=X # number of the technical replicate set e.g. 1
#	REP_LABEL=${LABEL}_rep$REP_NUM
#	INPUTS_TECH_REPS=('bash' 'array' 'of' 'bams' 'from' 'replicates') # BAM files you want combined as technical replicates
#   example bash array - INPUTS_TECH_REPS=('INPUT=A.L1.bam' 'INPUT=A.L2.bam' 'INPUT=A.L3.bam')
#	java -Xmx8G -Djava.io.tmpdir=temp/ -jar $PICARD MergeSamFiles $INPUTS_TECH_REPS OUTPUT=$TMP_DIR/$REP_LABEL.bam USE_THREADING=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT

echo "### Step 4: Mark duplicates"
#java -Xmx30G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=temp/ -jar $PICARD MarkDuplicates INPUT=$PAIR_DIR/$SRA.bam OUTPUT=$REP_DIR/$REP_LABEL.bam METRICS_FILE=$REP_DIR/metrics.$REP_LABEL.txt TMP_DIR=$TMP_DIR ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

java -Xmx30G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=tmp/ -jar $PICARD MarkDuplicates INPUT=$PAIR_DIR/$SRA.bam OUTPUT=$REP_DIR/$REP_LABEL.bam METRICS_FILE=$REP_DIR/metrics.$REP_LABEL.txt TMP_DIR=$TMP_DIR ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

$SAMTOOLS index $REP_DIR/$REP_LABEL.bam

perl $STATS $REP_DIR/$REP_LABEL.bam > $REP_DIR/$REP_LABEL.bam.stats

echo "Finished Mapping Pipeline through Duplicate Removal"

#########################################################################################################################################
###                                       How to Accommodate Biological Replicates                                                    ###
### This pipeline is currently built for processing a single sample with one read1 and read2 FASTQ file.                              ###
### Biological replicates (eg. multiple libraries made from the same sample) should be merged before proceeding with subsequent steps.###
### The code below is an example of how to merge biological replicates.                                                               ###
#########################################################################################################################################
#
#	INPUTS_BIOLOGICAL_REPS=('bash' 'array' 'of' 'bams' 'from' 'replicates') # BAM files you want combined as biological replicates
#   example bash array - INPUTS_BIOLOGICAL_REPS=('INPUT=A_rep1.bam' 'INPUT=A_rep2.bam' 'INPUT=A_rep3.bam')
#
#	java -Xmx8G -Djava.io.tmpdir=temp/ -jar $PICARD MergeSamFiles $INPUTS_BIOLOGICAL_REPS OUTPUT=$MERGE_DIR/$LABEL.bam USE_THREADING=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT
#
#	$SAMTOOLS index $MERGE_DIR/$LABEL.bam

# perl $STATS $MERGE_DIR/$LABEL.bam > $MERGE_DIR/$LABEL.bam.stats

# echo "Finished Mapping Pipeline through merging Biological Replicates"

RE='"GATC,GANTC,CTNAG,TTAA"'
HIC_BAM="$TMP_DIR3/$SRA.bam"

source /opt/anaconda3/etc/profile.d/conda.sh
conda activate haphic
### Run HapHiC pipeline
mkdir ${NAME}_HapHIC
cd ${NAME}_HapHIC
#srun -c 100 /data2/castri/HapHiC/haphic pipeline ${ASSEMBLY} ${HIC_BAM} ${CHR} --correct_nrounds 2 --remove_allelic_links ${PLOIDY} --threads 100 --RE $RE --prefix ${NAME}_Haphic

srun -c 100 ~/HapHiC/haphic pipeline ${ASSEMBLY} ${HIC_BAM} ${CHR} --correct_nrounds 2 --remove_allelic_links ${PLOIDY} --threads 100 --RE $RE --prefix ${NAME}_Haphic

### Create files for Juicebox
cd ./04.build 
bash juicebox.sh

cd ../..

conda deactivate
done
