export WORKSPACE=${WORKSPACE}
export REGWORK=/discover/nobackup/modele
export MODELROOT=/discover/nobackup/modele/devel/master
export MODELEBASELINE=/discover/nobackup/modele/modelE_baseline
export GCMSEARCHPATH=/discover/nobackup/projects/giss/prod_input_files
/discover/nobackup/modele/devel/master/exec/testing/modelEunitTests.sh
mail -s "modelE unit tests results" ccruz@nccs.nasa.gov < /home/modele/exec/modelEut-email.log
mail -s "modelE unit tests results" tclune@nccs.nasa.gov < /home/modele/exec/modelEut-email.log
