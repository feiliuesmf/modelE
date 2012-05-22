export WORKSPACE=${WORKSPACE}

# The following variables are used by mainRegTests.sh:
export REGWORK=/discover/nobackup/modele
export MODELEBASELINE=/discover/nobackup/ccruz/modelE_baseline
export GCMSEARCHPATH=/discover/nobackup/projects/giss/prod_input_files
export MODELROOT=/discover/nobackup/modele/devel/master.mock

# Clean up "signal" file in workspace
rm -f ${WORKSPACE}/.success

# This is where modelE testing scripts reside 
cd $MODELROOT/exec/testing

CONFIG=regTest
# Create a configuration file for this run:
rm -f $CONFIG.cfg
cat << EOF > $CONFIG.cfg
# decks arrays lists the rundecks to test:
@decks = ('EM20', 'E1oM20','E4TcadF40','E4F40','E4arobio_h4c','E4arobio_g6c');
# Git branch:
\$gitbranch = 'master';
# comps arrays specifies the compilers to use:
@comps = ('intel', 'gfortran');
# This specifies the test level. Other options are AGGRESSIVE and INSANE
\$level = 'GENTLE';
# Do we want to clean up the scratch space?
\$doCleanScratch = 'YES';
EOF

# This tells the scripts that this is a mock run
export MOCKMODELE=1

# Execute main test script 
mainRegTests.sh $CONFIG

chmod g+rw $MODELROOT/exec/testing/*.o*

# If everything is OK then scripts will write a "signal" file to workspace
if [ -e ${WORKSPACE}/.success ]; then 
   exit 0 
else 
   exit 1
fi
