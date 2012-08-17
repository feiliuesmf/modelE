#!/usr/local/bin/bash

# Script to run modelE unit tests on Linux (DISCOVER)

# -------------------------------------------------------------------
# FUNCTIONS
# -------------------------------------------------------------------

# -------------------------------------------------------------------
watchJob()
# -------------------------------------------------------------------
{
# Monitor job
# Input arguments: $1=job id
    jobID=$1

    maxWait=900
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

  local compiler=$1
  local jobScript=$2
  local testLog=$3

  local deck=E4TcadF40

# CC: This condition is only used to test the script.
  local mpi=YES
  if [ "$mpi" == "NO" ]; then
    pfunitSuffix=".serial"
  fi

# CREATE JOB SCRIPT

  cat << EOF > $jobScript
#!/bin/bash
#PBS -N modelEut
#PBS -l select=1:ncpus=12
#PBS -l walltime=0:15:00
#PBS -W group_list=s1001
#PBS -j oe
#PBS -V

# set up the modeling environment
. /usr/share/modules/init/bash
module purge
EOF

   if [ "$compiler" == "intel" ]; then

   cat << EOF >> $jobScript
module load comp/intel-12.1.0.233 mpi/impi-3.2.2.006
EOF

   else

   cat << EOF >> $jobScript
module load other/comp/gcc-4.7-20120331 other/mpi/mvapich2-1.8a2/gcc-4.7-20120331
EOF

   fi

   cat << EOF >> $jobScript

export PFUNIT=/discover/nobackup/modele/libs/pFUnit.${compiler}${pfunitSuffix}
export MODELERC=$REGSCRATCH/${compiler}/modelErc.${compiler}

# CC: This condition is only used to test the script.
if [ ! -e "\$MODELERC" ]; then
  export MODELERC=$HOME/.modelErc.${compiler}
fi

cd $REGSCRATCH
git clone $MODELROOT ${deck}.${compiler}

cd $REGSCRATCH/${deck}.${compiler}/decks
make rundeck RUN=$deck RUNSRC=$deck
EOF

   if [ "$compiler" == "intel" ]; then

   cat << EOF >> $jobScript
make -j gcm RUN=$deck EXTRA_FFLAGS="-O0 -g -traceback" MPI=$mpi
EOF

   else

   cat << EOF >> $jobScript
make -j gcm RUN=$deck EXTRA_FFLAGS="-O0 -g -fbacktrace" MPI=$mpi
make tests RUN=$deck MPI=$mpi >> $testLog 2>&1
EOF

   fi

   cat << EOF >> $jobScript
make tests RUN=$deck MPI=$mpi >> $testLog 2>&1
unset PFUNIT
EOF
  chmod +x $jobScript

# SUBMIT JOB SCRIPT

   echo "RESULTS [$compiler]:" >> $toEmail
   echo ""  >> $toEmail
   if [ "$mpi" == "YES" ]; then
     jobID=`qsub $jobScript`
     jobID=`echo $jobID | sed 's/.[a-z]*$//g'`
     if [ -z "$jobID" ]; then
       echo "There was a queue submission problem" >> $toEmail
       echo ""  >> $toEmail
       return
     fi

     watchJob $jobID jobRan

      if [ $jobRan -eq 0 ]; then
        echo "The PBS $jobID did not complete on time. Check $testLog" >> $toEmail
        echo ""  >> $toEmail
        return
      fi  
    else
      ./$jobScript
      wait
    fi

# PARSE FOR ERRORS

   local linesToShow=2
   local results=`grep 'run,.* failed' $testLog`
   local nbrOfFailed=`echo $results | sed 's/^[0-9]*\s*run, //g' | sed 's/\s*failed\s*.*$//g'`

   local anyError=`grep '[tests] Error' $testLog`
   if [ -n "$anyError" ]; then
     local errMsg="Error detected while compiling/running unit tests"
     linesToShow=10
   fi

   if [ "$anyError" != "" ]; then
     echo "$errMsg"  >> $toEmail
     echo "SUMMARY:"  >> $toEmail
     echo "..."  >> $toEmail
     tail -$linesToShow $testLog >> $toEmail
     echo ""  >> $toEmail

   else

     tail -$linesToShow $testLog >> $toEmail
     echo ""  >> $toEmail

  fi 

# DELETE JOB SCRIPT
  rm -f $jobScript

}

# ---------------------
# MAIN
# ---------------------

ROOT=$MODELROOT/exec/testing/
cd $ROOT
toEmail="$CONFIG.unit"
rm -f $toEmail
compilers=(intel gfortran)
for compiler in "${compilers[@]}"; do 
  job=modelE.${compiler}.j
  log=${ROOT}${compiler}".log"
  submitJob "$compiler" "$job" "$log"
  rm -f $job $log
done 

exit 0
