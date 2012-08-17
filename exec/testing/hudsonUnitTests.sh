#!/bin/bash
export WORKSPACE=${WORKSPACE}
# Path of scratch space
export REGWORK=/discover/nobackup/ccruz
# This is where the reference restarts are kept
export MODELEBASELINE=/discover/nobackup/modele/modelE_baseline
# modelE data is here
export GCMSEARCHPATH=/discover/nobackup/projects/giss/prod_input_files
# Path of up-to-date modelE repository
export MODELROOT=/discover/nobackup/modele/devel/modelE.clones/master

$MODELROOT/exec/testing/modelEunitTests.sh
mail -s "modelE unit tests results" ccruz@nccs.nasa.gov < /home/modele/exec/modelEut-email.log
mail -s "modelE unit tests results" tclune@nccs.nasa.gov < /home/modele/exec/modelEut-email.log
