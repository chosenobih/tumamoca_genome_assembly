### First assembly assessment script. 
### Meant to be run under run_assembly_assessment.sh 

ASSEMBLY=$1
BUSCO=$2
ASSEMBLY_PATH=$3
DATE=$(date +'%d.%m.%Y')

mkdir Assembly_Assessment_${ASSEMBLY}

### Run assembly stats
assembly-stats ${ASSEMBLY_PATH}/${ASSEMBLY}.fa > ./Assembly_Assessment_${ASSEMBLY}/${ASSEMBLY}_assembly_stats.txt

### Run BUSCO
srun -c 64 busco -i ${ASSEMBLY_PATH}/${ASSEMBLY}.fa -l ${BUSCO} -c 64 -o ./Assembly_Assessment_${ASSEMBLY}/${ASSEMBLY} -m genome

