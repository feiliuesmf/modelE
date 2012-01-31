#!/usr/local/bin/bash

# This is the top level script to run modelE regression tests.
# It is invoked by a cron job and runs the perl scripts under modelE/exec.

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

   cd $NOBACKUP/devel/master/exec
   echo "Execute regressionTests.pl..."
   /usr/bin/perl regressionTests.pl > nohup.out 2>&1
   wait
   jobID=`qsub $NOBACKUP/devel/master/exec/diffreport.j`
   jobID=`echo $jobID | sed 's/.[a-z]*$//g'`
   watch_job $jobID

mail -s "discover results" giss-modelE-regression@lists.nasa.gov < $NOBACKUP/devel/master/exec/DiffReport

echo "Done".
