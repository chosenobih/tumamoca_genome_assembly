### Third script of assembly assessment pipeline
### Meant to be run under run_assembly_assessment.sh

ASSEMBLY=$1
READS=$2
DATE=$(date +'%d.%m.%Y')
ASSEMBLY_PATH=$3

### Plot Coverage

### Map Reads
bwa index ${ASSEMBLY_PATH}/${ASSEMBLY}.fa 
srun -c 128 minimap2 -I 128G -ax map-pb -t 128 -c --secondary=no ${ASSEMBLY_PATH}/${ASSEMBLY}.fa ${READS} | samtools view -bh > ${ASSEMBLY_PATH}/HiFi2_${ASSEMBLY}.bam
samtools sort ${ASSEMBLY_PATH}/HiFi2_${ASSEMBLY}.bam -o ${ASSEMBLY_PATH}/HiFi2_${ASSEMBLY}_sorted.bam
samtools index ${ASSEMBLY_PATH}/HiFi2_${ASSEMBLY}_sorted.bam

### Create Dictionary
picard CreateSequenceDictionary -R ${ASSEMBLY_PATH}/${ASSEMBLY}.fa

module load java/23

### Plot: may need to edit contig length
java -jar /data2/castri/tools/jvarkit/jvarkit.jar wgscoverageplotter --dimension 5000x1000 -C -1 --clip -R ${ASSEMBLY_PATH}/${ASSEMBLY}.fa ${ASSEMBLY_PATH}/HiFi2_${ASSEMBLY}_sorted.bam --percentile median --min-contig-length 1000000 >  ./Assembly_Assessment_${ASSEMBLY}/${ASSEMBLY}.svg
