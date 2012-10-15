export WORKSPACE=${WORKSPACE}

# The following variables are used by mainRegTests.sh:
export RUN_UNIT_TESTS=NO
export RUN_TESTS=YES
export CREATE_DIFF=NO
export REGWORK=/discover/nobackup/modele
export MODELEBASELINE=/discover/nobackup/modele/modelE_baseline
export GCMSEARCHPATH=/discover/nobackup/projects/giss/prod_input_files
export MODELROOT=/discover/nobackup/modele/devel/master

# Clean up "signal" file in workspace
rm -f ${WORKSPACE}/.success

# This is where modelE testing scripts reside 
cd $MODELROOT/exec/testing

CONFIG=regTest
# Create a configuration file for this run:
rm -f $CONFIG.cfg
cat << EOF > $CONFIG.cfg
# The configurations hash consists of a rundeck and its associated options
# Format is as follows: rundeck => [options], where rundeck is a rundeck that
# belongs to the rundeck suite and [options] are entered in the following
# order:
#   [options] = mode, <durations array>, debug
#   where
#     mode         : S(SERIAL), M(MPI), B(BOTH)
#     durations    : If it equals to 0 (default) then it runs for 
#                    1hr+1dy+restart (used for testing)
#                    Else entry it runs for the specified time 
#                    (in model time-steps, e.g. 96=2dy, 1440=1mo)
#    debug option  : Y or N (allows for compilation with debug flags)

%configurations = ( 'EM20'   => ['B',0,'N'],
                    'E1oM20' => ['B',0,'N'] );

#                    'E4F40'     => ['B',0,'N']
#                    'E4TcadF40' => ['M',0,'N']
#                    'E4arobio_g6c' => ['M',0,'N']
#                    'E4arobio_h4c' => ['M',0,'N']
#                    'E4C90L40' => ['M',0,'N']
#                    'SCMSGPCONT' => ['S',0,'N']
#                    'E_AR5_CADI' => ['M',0,'N']
#                    'nonProduction_E4TcadC12' => ['B',1440,'N']
#                    'nonProduction_E_AR5_C12' => ['B',1440,'N']
#                    'E1oM20' => ['M',0,'N']

# Specify the compilers to use: intel, gfortran, nag
@compilers = ('intel','gfortran');

# Git branch to clone and test:
\$gitbranch = 'master';

# This specifies the test level. Other options are AGGRESSIVE and INSANE. 
# Since modelE domain decomposition is along latitudinal direction the
# different options vary the number of PEs as a function of NLATS.
\$level = 'GENTLE';

# Do we want to clean up the scratch space? Generally yes.
\$doCleanScratch = 'YES';
EOF

# Execute main test script 
mainRegTests.sh $CONFIG

# Make PBS output files readable to group
chmod g+rw *.o*

# If everything is OK then scripts will write a "signal" file to workspace
if [ -e ${WORKSPACE}/.success ]; then 
   exit 0 
else 
   exit 1
fi
