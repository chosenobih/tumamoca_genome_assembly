NAME=$1
HiFi=$2
HIC_R1=$3
HIC_R2=$4

echo "Processing ${ASSEMBLY}"
#### Edit ploidy here if needed 
#l_value="3"

echo $HiFi

###Activate conda env for packages
#conda activate rice_pang		#activate conda env before running the script

###run KMC:
mkdir ${NAME} ${NAME}_kmc
srun -c 120 kmc -k21 -m100 -t120 ${HiFi} ${NAME} ${NAME}_kmc
kmc_tools -t120 transform ${NAME} histogram ${NAME}_hist.txt

### run hifiasm
echo "Starting assembly steps"
mkdir ${NAME}_hifiasm
srun -c 200 hifiasm -t 200 -o ${NAME}_hifiasm --dual-scaf --telo-m TTTAGGG --primary -l3 --h1 "$HIC_R1" --h2 "$HIC_R2" --write-ec ${HiFi} 2> ./${NAME}_hifiasm.stdout
mv ${NAME}_hifiasm.* ./${NAME}_hifiasm

### Fix output gfas to fas
awk '/^S/{print ">"$2;print $3}' ./${NAME}_hifiasm/${NAME}_hifiasm.hic.p_ctg.gfa > ./${NAME}_hifiasm/${NAME}_hifiasm.hic.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' ./${NAME}_hifiasm/${NAME}_hifiasm.hic.a_ctg.gfa > ./${NAME}_hifiasm/${NAME}_hifiasm.hic.a_ctg.fa
awk '/^S/{print ">"$2;print $3}' ./${NAME}_hifiasm/${NAME}_hifiasm.hic.hap1.p_ctg.gfa  > ./${NAME}_hifiasm/${NAME}_hifiasm.hic.hap1.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' ./${NAME}_hifiasm/${NAME}_hifiasm.hic.hap2.p_ctg.gfa  > ./${NAME}_hifiasm/${NAME}_hifiasm.hic.hap2.p_ctg.fa
#cat ./${NAME}_hifiasm/${NAME}_hifiasm.p_ctg.fa ./${NAME}_hifiasm/${NAME}_hifiasm.a_ctg.fa > ./${NAME}_hifiasm/${NAME}_hifiasm.combined_ctg.fa

### Map reads to prepare for polishing primary and combined haplotype assembly. Hash out if only one is needed.
srun -c 120 minimap2 -x map-pb -I 160G -a -t 120 -c --secondary=no ./${NAME}_hifiasm/${NAME}_hifiasm.hic.p_ctg.fa ${HiFi} > ./${NAME}_hifiasm/${NAME}_HiFi_p.sam

### run racon polishing
mkdir ${NAME}_racon
srun -c 120 racon -t 120 ${HiFi} ./${NAME}_hifiasm/${NAME}_HiFi_p.sam ./${NAME}_hifiasm/${NAME}_hifiasm.hic.p_ctg.fa > ./${NAME}_racon/${NAME}_hifiasm.p_racon.fa

### Now assembly assessment
ASSEMBLY="${NAME}_hifiasm.p_racon"
ASSEMBLY_PATH="./${NAME}_racon"
BUSCO="embryophyta_odb10"

### Run assembly_stats and BUSCO
sh ./scripts/assembly_assessment_1.sh ${ASSEMBLY} ${BUSCO} ${ASSEMBLY_PATH}

### Run Merqury
## May need to edit range for plot
sh ./scripts/assembly_assessment_2.sh ${ASSEMBLY} ${HiFi} ${ASSEMBLY_PATH}


#### Map reads and create coverage plot of assembly. 
### May need to edit plotting length requirements for in assembly_assessment_3.sh if outside of range
sh ./scripts/assembly_assessment_3.sh ${ASSEMBLY} ${HiFi} ${ASSEMBLY_PATH}

### Files housekeeping
mv ${NAME}* ${NAME}/
mv log* ${NAME}/
mv Assembly_Assessment_${NAME}* ${NAME}

echo "####################################################################"
echo "Pipeline run done"
date
echo "####################################################################"
