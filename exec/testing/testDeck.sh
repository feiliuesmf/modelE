#!/bin/bash

# This script is invoked from HUDSON. It takes in an argument:
# the rundeck configuration file name (all lowercase, no extension)
# It sets ENV variables used by the regression scripts and finally
# it runs the main regression test script.

deck=$1

echo "Exported HUDSON ENVS:"
echo "WORKSPACE="$WORKSPACE

# The following variables are used by mainRegTests.sh:
export REGWORK=/discover/nobackup/modele
export MODELEBASELINE=/discover/nobackup/modele/modelE_baseline
export GCMSEARCHPATH=/discover/nobackup/projects/giss/prod_input_files
export MODELROOT=/discover/nobackup/modele/devel/master

cd $MODELROOT/exec/testing

mainRegTests.sh "$deck"

# Change permissions of PBS output files
chmod g+rw *.o*
