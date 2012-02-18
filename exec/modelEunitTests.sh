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
       mail -s "$subject" giss-modelE-regression@lists.nasa.gov < "$toEmail"
    else
       subject="modelE unit tests: failed"
       echo -e "$subject" | mail -s "$subject" giss-modelE-regression@lists.nasa.gov < "$testLog"
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

    cat << EOF > $work_path/modelEut.j
#!/bin/bash
#PBS -N modelEut
#PBS -l select=1
#PBS -l walltime=1:00:00
#PBS -W group_list=k3002
#PBS -j oe
#PBS -V

# set up the modeling environment
. /usr/share/modules/init/bash
module purge
module load comp/intel-11.1.072 mpi/impi-3.2.2.006

export PFUNIT=/discover/nobackup/modele/libs/pFUnitIntel
export MODELERC=/discover/nobackup/modele/regression_scratch/intel/modelErc.intel

cd /discover/nobackup/modele/regression_scratch
git clone /discover/nobackup/modele/regression_scratch/master pFUnit

cd /discover/nobackup/modele/regression_scratch/pFUnit/decks
make rundeck RUN=unitTests RUNSRC=E4TcadF40
make -j gcm RUN=unitTests EXTRA_FFLAGS="-O0 -g -traceback" MPI=YES
cd /discover/nobackup/modele/regression_scratch/pFUnit/tests
make tests > $testLog 2>&1
wait
tail -2 $testLog >> $toEmail
EOF

    jobID=`qsub $work_path/modelEut.j`
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
}

# ---------------------
# MAIN
# ---------------------

work_path=/home/modele/exec
toEmail=$work_path"/modelEut-email.log"
testLog=$work_path"/modelEut-tests.log"
cd $work_path
rm -f $toEmail $testLog 

submitJob
sendEmailReport

exit 0
