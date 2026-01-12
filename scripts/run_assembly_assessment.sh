### Run assembly assessment statistics
#### Requirements: BUSCO (check env name and edit before running script), R library with Merqury dependencies installed, assembly-stats, wgsplotter
#### To run- Edit paths/names below then: source run_assembly_pipelin.sh
#### Associated assembly assessment scripts must be in the same folder

########## Change file names here before running ##########
ASSEMBLY=$1
READS=$2
ASSEMBLY_PATH=$3
BUSCO="embryophyta_odb10"

#### Run assembly_stats and BUSCO
sh ./scripts/assembly_assessment_1.sh ${ASSEMBLY} ${BUSCO} ${ASSEMBLY_PATH}

#### Run Merqury
sh ./scripts/assembly_assessment_2.sh ${ASSEMBLY} ${READS} ${ASSEMBLY_PATH}

#### Map reads and create coverage plot of assembly. 
#### May need to edit plotting length requirements for in assembly_assessment_3.sh if outside of range
sh ./scripts/assembly_assessment_3.sh ${ASSEMBLY} ${READS} ${ASSEMBLY_PATH}
