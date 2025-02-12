#!/bin/bash

# Installation
### build singularity image
#IMAGE="immcantation_suite-4.5.0.sif"
#singularity build $IMAGE docker://immcantation/suite:4.5.0
### re-build with apptainer
#apptainer build immcantation_suite-4.5.0.sif docker://immcantation/suite:4.5.0
### create a symbolic link of this sif file to which ever directory you want to run the apptainer command from
### to run the immantation image, the sif file must be in your current directory

### TEMPLATE: run in terminal (running this in screen)
#           - using singularity/apptainer, the paths you specificy are to locations on the server 
#             not within the immcantation image
#           - use full file and directory PATHs, no relative PATHs
# apptainer run immcantation_suite-4.5.0.sif [PATH]/impact-analysis/code/0_vdj_alignments/2_run_immcantation.bash [PATH]/input_immcantation_one_patient [PATH]/results_immcantation_one_patient


# Arguments
# - Input formatted sample name
DATA_DIR=$1
OUT_DIR=$2

READS_BCR=$DATA_DIR\/input_bcr_filtered_contig.fasta
READS_TCR=$DATA_DIR\/input_tcr_filtered_contig.fasta
ANNOTATIONS_BCR=$DATA_DIR\/input_bcr_filtered_contig_annotations.csv
ANNOTATIONS_TCR=$DATA_DIR\/input_tcr_filtered_contig_annotations.csv
OUT_DIR_BCR=$OUT_DIR\/BCR_CSFPB
OUT_DIR_TCR=$OUT_DIR\/TCR_CSFPB
MODEL=aa
NPROC=20

BCR_DIST=0.15  # tested on full cohort

# >>> BCR
changeo-10x -s $READS_BCR -a $ANNOTATIONS_BCR -x $BCR_DIST -n BCR_CSFPB -o $OUT_DIR_BCR -p $NPROC -t ig -m $MODEL -f airr

# >>> TCR
#changeo-10x -s $READS_TCR -a $ANNOTATIONS_TCR -x 0.00 -n TCR_CSFPB -o $OUT_DIR_TCR -p $NPROC -t tr -m $MODEL -f airr


# >>> process TCR Alpha chain results
#cd $OUT_DIR_TCR
#DefineClones.py -d TCR_CSFPB_light_productive-T_sc-pass.tsv --act first --model $MODEL --dist 0.00
#CreateGermlines.py -d TCR_CSFPB_light_productive-T_sc-pass_clone-pass.tsv -o TCR_CSFPB_light_germ-pass.tsv -g dmask --cloned -r /usr/local/share/germlines/imgt/human/vdj/imgt_human_TRA*

# >>> process BCR light chain results
cd $OUT_DIR_BCR
DefineClones.py -d BCR_CSFPB_light_productive-T_sc-pass.tsv --act first --model $MODEL --dist $BCR_DIST
CreateGermlines.py -d BCR_CSFPB_light_productive-T_sc-pass_clone-pass.tsv -o BCR_CSFPB_light_germ-pass.tsv -g dmask --cloned -r /usr/local/share/germlines/imgt/human/vdj/imgt_human_IG*

