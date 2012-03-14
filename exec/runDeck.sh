#!/bin/bash -f

# runDeck: This script runs modelE rundecks
# History: C.Cruz - created 1/2011 

# -------------------------------------------------------------------
# MAIN
# -------------------------------------------------------------------

   arch=`uname -s`
   node=`uname -n`

# If running on DISCOVER restrict execution of this script to be under PBS
# environment INTERACTIVE or BATCH

   if [[ "$node" =~ discover || "$node" =~ dali ]]; then
      [ -z "$PBS_ENVIRONMENT" ] && 
      echo " ### Need to run under PBS interactive or batch job." && 
      exit 1;
   fi
   
   if [[ "$node" =~ borg ]]; then
      pbsType=`echo $PBS_ENVIRONMENT`
   fi

# runDeckFunctions.sh must exists
   command -v runDeckFunctions.sh &>/dev/null || 
   { 
      echo " ~~~ $0 cannot run: runDeckFunctions.sh is missing." >&2
      exit 1 
   }
   . runDeckFunctions.sh

# Number of command line arguments
   numargs=$#
   if [ $numargs -lt 1 ]; then
      echo " ### Please specify git repository"; exit 1
   fi

# define defaults
   defaults

# get options from command line arguments
   OIFS=$IFS; IFS=$(echo -e "\n"); getopt names; IFS=$OIFS
   names=("$@")

# parse command line arguments
   parseOptions

# Initialize settings
   initEnvironment
   
# Build/run the deck
   buildAndRun

# compare run restarts
   compareRuns
   if [ "$restartRegression" == "YES" ]; then
      compareRegressionRuns
   fi

# clean up work space
   finalize 0
