### Identify and extract organelles from assembly then incorporate into final assembly 
### Requires: oatk and appropriate organelle database (in working directory),
### Scaffolded FASTA file  and name of run
### To run sh add_organelles.sh <scaffolded fasta file> <output files prefix name>
ASSEMBLY=$1
NAME=$2
Chr=14

######## Edit info below #########
OUTPUT="${NAME}"
MITO="./angiosperm_mito.fam"
CP="./angiosperm_pltd.fam"


samtools faidx ${ASSEMBLY}

sort -nk 1 ${ASSEMBLY}.fai | cut -f 1 | head -n ${Chr} > primary_scaffolds.list
scaffold_count=`wc -l ${ASSEMBLY}.fai |  awk '{ print $1 }'`
unanchored_count=`echo $((${scaffold_count}-${Chr}))`
sort -nk 1 ${ASSEMBLY}.fai | cut -f 1 | tail -n ${unanchored_count} > unanchored_scaffolds.list

seqtk subseq ${ASSEMBLY} primary_scaffolds.list > ./${NAME}_primary.fa
seqtk subseq ${ASSEMBLY} unanchored_scaffolds.list > ./${NAME}_unanchored.fa

CHR_ASSEMBLY="./${NAME}_primary.fa"
UNANCHORED="./${NAME}_unanchored.fa"

oatk -k 41 -c 15 -t 8 -m ${MITO} -p ${CP} -o ${OUTPUT} ${UNANCHORED}

######  If there is an issue with the pipeline, these are subscripts that can be used to initate the later stages
#hmm_annotation -t 8 -o ${OUTPUT}.annot_mito.txt ${MITO} ${OUTPUT}.utg.final.gfa
#hmm_annotation -t 8 -o ${OUTPUT}.annot_pltd.txt ${CP} ${OUTPUT}utg.final.gfa
#pathfinder -m ${OUTPUT}.annot_mito.txt -p ${OUTPUT}.annot_pltd.txt -o ${OUTPUT} ${OUTPUT}.utg.final.gfa

cat ${CHR_ASSEMBLY} ${OUTPUT}.mito.ctg.fasta  ${OUTPUT}.pltd.ctg.fasta > ${OUTPUT}_primary_organelles.fasta

## Map back unanchored scaffolds to organelles to extract mt and pltd scaffolds
minimap2 -I 128G -ax map-pb -t 20 -c --secondary=no ${OUTPUT}_primary_organelles.fasta ${UNANCHORED} > ${OUTPUT}_unanchored2primary_organelle.sam

### Get names of organelles from oatk
grep ">" ${OUTPUT}.mito.ctg.fasta | sed 's/>//g' | cut -f 1 > mito.list
grep ">" ${OUTPUT}.pltd.ctg.fasta | sed 's/>//g' | cut -f 1 > pltd.list
cat mito.list pltd.list > organelles.list

### Prepare bam file
samtools sort ${OUTPUT}_unanchored2primary_organelle.sam -o ${OUTPUT}_unanchored2primary_organelle.bam
samtools index ${OUTPUT}_unanchored2primary_organelle.bam

### Extract scaffolds mapping to organelles
samtools view -q 20 ${OUTPUT}_unanchored2primary_organelle.bam -R organelles.list > ${OUTPUT}_organelles.sam
### Get list of scaffolds mapped to organelles
cut -f 1 ${OUTPUT}_organelles.sam > organelle_scaffolds.list
### Get list of all initial assembly scaffolds
grep ">" ${ASSEMBLY} | sed 's/>//g' | cut -f 1 | sort > all_scaffolds.list
#### Remove organelle scaffolds from all assembly scaffolds so they are not in assembly twice
sort organelle_scaffolds.list -o organelle_scaffolds_sorted.list
comm -3 all_scaffolds.list organelle_scaffolds_sorted.list > scaffolds_no_organelles.list

### Create assembly FASTA without organelles or the scaffolds aligned to them
seqtk subseq ${ASSEMBLY} scaffolds_no_organelles.list > ${ASSEMBLY}_no_organelles.fa

### Combine assembly with no organelles with the oatk FASTA files of organelles for final assembly
### Stupid name sorry, rename accordingly
cat ${ASSEMBLY} ${OUTPUT}.mito.ctg.fasta ${OUTPUT}.pltd.ctg.fasta > ${OUTPUT}_final_with_organelles.fa
