#!/bin/bash

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

# runModelEfuncs.sh must exists
   command -v runModelEfuncs.sh &>/dev/null || 
   { 
      echo " ~~~ $0 cannot run: runModelEfuncs.sh is missing." >&2
      exit 1 
   }
   . runModelEfuncs.sh

   if [ "$#" -eq 1 ]; then
      if [[ "$1" =~ "-h" ]]; then
         usage; exit 0
         exit 1
      fi      
   fi      

# define defaults
   defaults

# get options from command line arguments
   OIFS=$IFS; IFS=$(echo -e "\n"); getopt names; IFS=$OIFS
   names=("$@")

# parse command line arguments
   parseOptions

# Required ENVIROMENT variables
   if [ -z "$MODELEGIT" ]; then
      echo " ### Please specify modelE git repository: e.g. export MODELEGIT=<path>"; exit 1
   else 
      gitRepository=$MODELEGIT
   fi

# Initialize settings
   initEnvironment

# Build/run the deck
   buildAndRun
   if [[ "$check" == "NO"  && "$baseline" == "YES" ]]; then finalize 0; fi

# compare run restarts with baseline
   compareRuns

# clean up work space
   finalize 0
