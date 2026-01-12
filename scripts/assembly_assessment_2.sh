### Second script of assembly assessment pipeline
### Meant to be run under run_assembly_assessment.sh

ASSEMBLY=$1
READS=$2
ASSEMBLY_PATH=$3
DATE=$(date +'%d.%m.%Y')

### Run Merqury (it can be a diva)

srun -c 128 meryl count threads=128 k=21 ${READS} output ${READS}.meryl
ln -s $MERQURY/merqury.sh

### It has issues if the assembly file isn't in the working directory

cp ${ASSEMBLY_PATH}/${ASSEMBLY}.fa .
srun -c 128 $MERQURY/merqury.sh ${READS}.meryl ${ASSEMBLY}.fa ${ASSEMBLY}_Merqury

### Adjust axes with -m and -n of plots 
Rscript /opt/marbl/merqury-1.3/plot/plot_spectra_cn.R -f ${ASSEMBLY}_Merqury.${ASSEMBLY}.spectra-cn.hist -o ${ASSEMBLY}_Merqury.spectra-cn.hist  -n 11119921 -m 200 

mv ${ASSEMBLY}_Merqury.* ./Assembly_Assessment_${ASSEMBLY}
mv *.hist ./Assembly_Assessment_${ASSEMBLY}
mv *.filt ./Assembly_Assessment_${ASSEMBLY}
mv *.ploidy ./Assembly_Assessment_${ASSEMBLY}
mv *.wig ./Assembly_Assessment_${ASSEMBLY}
