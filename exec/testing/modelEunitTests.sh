#!/usr/local/bin/bash

# Script to run modelE unit tests on Linux (DISCOVER)
# It is invoked by a cron job

# -------------------------------------------------------------------
# FUNCTIONS
# -------------------------------------------------------------------

# -------------------------------------------------------------------
gatherForEmail()
# -------------------------------------------------------------------
{
    echo $1 >> $toEmail 2>&1
}

# -------------------------------------------------------------------
sendEmailReport()
# -------------------------------------------------------------------
{
# email info
    if [ "$FAIL" == "" ]; then
       subject="modelE unit tests: successful"
#       mail -s "$subject" giss-modelE-regression@lists.nasa.gov < "$toEmail"
    else
       subject="modelE unit tests: failed"
#       echo -e "$subject" | mail -s "$subject" giss-modelE-regression@lists.nasa.gov < "$testLog"
    fi
}

# -------------------------------------------------------------------
failScript()
# -------------------------------------------------------------------
{
    FAIL=YES
    echo $(date +%H:%M:%S) ERROR: $1 >> $toEmail 2>&1
    sendEmailReport
    exit 1
}

# -------------------------------------------------------------------
watch_job()
# -------------------------------------------------------------------
{
# Monitor job
# Input arguments: $1=job id
    jobID=$1

    maxWait=3600
    seconds=0
    jobSuccess=0
    while [ $seconds -lt $maxWait ];
      do
      qStatus=`qstat | grep $jobID | awk '{print $5}'`
      if [ -z "$qStatus" ]; then
          jobSuccess=1
          break
      fi
      sleep 10
      let seconds=$seconds+10
    done
    eval $2="$jobSuccess"
}
# -------------------------------------------------------------------
submitJob()
# -------------------------------------------------------------------
{

   compiler=$1
   jobScript=$work_path/modelE.${compiler}.j
   cat << EOF > $jobScript
#!/bin/bash
#PBS -N modelEut
#PBS -l select=1
#PBS -l walltime=1:00:00
#PBS -W group_list=s1001
#PBS -j oe
#PBS -V

# set up the modeling environment
. /usr/share/modules/init/bash
module purge
EOF

   if [ "$compiler" == "intel" ]; then

   cat << EOF >> $jobScript
module load comp/intel-12.1.0.233 mpi/impi-4.0.3.008
EOF

   else

   cat << EOF >> $jobScript
module load other/comp/gcc-4.7-20120331 other/mpi/mvapich2-1.8a2/gcc-4.7-20120331
EOF

   fi

   cat << EOF >> $jobScript

export PFUNIT=/discover/nobackup/modele/libs/pFUnit.${compiler}
export MODELERC=/discover/nobackup/modele/regression_scratch/${compiler}/modelErc.${compiler}

cd /discover/nobackup/modele/regression_scratch
git clone /discover/nobackup/modele/regression_scratch/master unitTests.${compiler}

cd /discover/nobackup/modele/regression_scratch/unitTests.${compiler}/decks
make rundeck RUN=unitTests RUNSRC=E4TcadF40
EOF

   if [ "$compiler" == "intel" ]; then

   cat << EOF >> $jobScript
make -j gcm RUN=unitTests EXTRA_FFLAGS="-O0 -g -traceback" MPI=YES
EOF

   else

   cat << EOF >> $jobScript
make -j gcm RUN=unitTests EXTRA_FFLAGS="-O0 -g" MPI=YES
EOF

   fi

   cat << EOF >> $jobScript
cd /discover/nobackup/modele/regression_scratch/unitTests.${compiler}/tests
EOF

   if [ "$compiler" == "intel" ]; then

   cat << EOF >> $jobScript
make tests >> $testLog 2>&1
EOF

   else

   cat << EOF >> $jobScript
make tests F90_VENDOR=GFortran >> $testLog 2>&1
EOF

   fi

   jobID=`qsub $jobScript`
   jobID=`echo $jobID | sed 's/.[a-z]*$//g'`
   if [ -z "$jobID" ]; then
      FAIL=YES
      gatherForEmail "   +++There was a queue submission problem"
      return
   fi

   watch_job $jobID jobRan

   if [ $jobRan -eq 0 ]; then
      failScript "The PBS $jobID did not complete on time."
   fi

   linesToShow=2
   results=`grep 'run,.* failed' $testLog`
   nbrOfFailed=`echo $results | sed 's/^[0-9]*\s*run, //g' | sed 's/\s*failed\s*.*$//g'`

   anyError=`grep 'make: ***' $testLog`
   if [ -n "$anyError" ]; then
     errMsg="Failure building unit tests"
     linesToShow=10
   fi

   echo "RESULTS [$compiler]:" >> $toEmail
   echo ""  >> $toEmail
   if [ -n "$anyError" ]; then

     echo "$errMsg"  >> $toEmail
     echo "SUMMARY:"  >> $toEmail
     echo "..."  >> $toEmail
     tail -$linesToShow $testLog >> $toEmail
     echo ""  >> $toEmail

   else

     tail -$linesToShow $testLog >> $toEmail
     echo ""  >> $toEmail

  fi #anyError

}

# ---------------------
# MAIN
# ---------------------

work_path=/home/modele/exec
toEmail=$work_path"/modelEut-email.log"
testLog=$work_path"/modelEut-tests.log"
cd $work_path
rm -f $toEmail $testLog 

submitJob "intel"
submitJob "gfortran"
sendEmailReport

exit 0
