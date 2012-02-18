#!/usr/local/bin/bash

# This is the top level script to run modelE regression tests.
# It is invoked by cron, Hudson, or user and runs the perl scripts under $TESTROOT/$MODELROOT/exec.
# TODO: Define $MODELROOT referenced above. Still actually hardcoded to devel/master/

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
   if  [ "$(need_default $TESTROOT)" == "true" ] ; then
        export TESTROOT=$NOBACKUP
   fi

}


watch_job()
{
# Monitor job
# Input arguments: $1=job id
   jobID=$1

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

   # User may just want the help page
   if [ "$1" = "--help" -o "$1" = "-h" -o ! -z "$1" ]
   then
        ./regTestsHelp.sh $0 $1; exit 0 ;
   fi

   get_defaults
   cd $TESTROOT/devel/master/exec
   echo "Execute regressionTests.pl..."
   /usr/bin/perl regressionTests.pl > nohup.out 2>&1
   wait
   if [ -z $MOCKMODELE ]; then
     jobID=`qsub $TESTROOT/devel/master/exec/diffreport.j`
     jobID=`echo $jobID | sed 's/.[a-z]*$//g'`
     watch_job $jobID
     mail -s "discover results" giss-modelE-regression@lists.nasa.gov < $TESTROOT/devel/master/exec/DiffReport
   else
     $TESTROOT/devel/master/exec/diffreport.j
     mail -s "mock modelE results" meandrew@nccs.nasa.gov < $TESTROOT/devel/master/exec/DiffReport
   fi

   echo "Done".
