### Quick assembly contamination check

### NCBI quick contamination screen
### Requires fcs with gxdb installed
### Get taxon ID from NCBI

# To run: sh contamination_screen.sh <assembly file> <output name> <NCBI taxon ID>

ASSEMBLY=$1
OUT_DIR=$2
TAX_ID=$3

module load python37
python3 /data2/service_center_analyses/rice_pangenomics/fcs-gx-0.5.0/fcs.py screen genome --fasta ${ASSEMBLY} --out-dir ${OUT_DIR} --gx-db /data2/service_center_analyses/rice_pangenomics/fcs-gx-0.5.0/gxdb --tax-id ${TAX_ID} --generate-logfile T
