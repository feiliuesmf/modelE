#!/bin/bash

# This is the top level script to run modelE regression tests.
# It is invoked by cron, Hudson, or user and runs the perl scripts under $MODELROOT/master/exec.

# Set up working variables
need_default(){
   if [ -z "$1" ]; then
        echo "true"
        exit 0
   else
        echo "false"
        exit 1
   fi
}
get_defaults()
{
   if  [ "$(need_default $MODELROOT)" == "true" ] ; then
	export MODELROOT=$NOBACKUP/devel/modelE.clones/simplex
   fi

   if  [ "$(need_default $REGWORK)" == "true" ] ; then
	export REGWORK=$NOBACKUP
   fi

   export REGSCRATCH=$REGWORK/regression_scratch
   # create if necessary
   if [ ! -d "$REGSCRATCH" ]; then
   	mkdir $REGSCRATCH
   fi

   export REGRESULTS=$REGWORK/regression_results
   # create if necessary
   if [ ! -d "$REGRESULTS" ]; then
   	mkdir $REGRESULTS
   fi

   if  [ "$(need_default $MODELEBASELINE)" == "true" ] ; then
	export MODELEBASELINE=$NOBACKUP/modelE_baseline
   fi

   if  [ "$(need_default $GCMSEARCHPATH)" == "true" ] ; then
	export GCMSEARCHPATH=/discover/nobackup/projects/giss/prod_input_files
   fi
	
}

watch_job()
{
# Monitor job
# Input arguments: $1=job id
   local jobID=$1

   maxWait=3600
   seconds=0
   done=0
   while [ $seconds -lt $maxWait ];
     do
     qStatus=`qstat | grep $jobID | awk '{print $5}'`
     if [ -z "$qStatus" ]; then
         done=1
         break
     fi
     sleep 30
     let seconds=$seconds+30
   done
}

# MAIN PROGRAM

   # User may just want the help page
   #if [ "$1" = "--help" -o "$1" = "-h" -o ! -z "$1" ]
   #if [[ "$1" == "--help" || "$1" == "-h" ]]
   #then
   #     ./regTestsHelp.sh $0 $1; exit 0 ;
   #fi

   if [ $# -ne 1 ]; then
      echo "Usage: `basename $0` {cfgFile}"
      exit 65
   fi
   export CFG_NAME=$1
   cfgFile=$1.cfg
   export CFG_FILE=$cfgFile

   get_defaults
   echo "Execute regressionTests.pl with "$cfgFile
   rm -f $1.out
   /usr/bin/perl regressionTests.pl $CFG_FILE > $1.out 2>&1
   # Not sure why I still need to do this:
   chmod g+rw $1.out
   chgrp s1001 $1.out LOG
 
   wait
   if [ -z $MOCKMODELE ]; then
     jobID=`qsub $MODELROOT/exec/testing/diffreport.j`
     jobID=`echo $jobID | sed 's/.[a-z]*$//g'`
     watch_job $jobID
#     mail -s "discover results" giss-modelE-regression@lists.nasa.gov < $MODELROOT/exec/testing/DiffReport
      cp ${CFG_NAME}.diff $WORKSPACE
   else
     $MODELROOT/exec/testing/diffreport.j 
     mail -s "mock modelE results" ccruz@nccs.nasa.gov < $MODELROOT/exec/testing/${CFG_NAME}.diff
   fi

   echo "Done".
