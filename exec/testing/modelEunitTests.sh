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
  local mpi=$4

  local deck=E4TcadF40

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
module load comp/intel-13.0.0.079 mpi/impi-3.2.2.006
EOF
   
  else

    cat << EOF >> $jobScript
module load other/comp/gcc-4.7-20120331 other/mpi/mvapich2-1.8a2/gcc-4.7-20120331
EOF

  fi

  cat << EOF >> $jobScript

export PFUNIT=/discover/nobackup/modele/libs/pFUnit.${compiler}${pfunitSuffix}
export MODELERC=$REGSCRATCH/${compiler}/modelErc.${compiler}

cd $REGSCRATCH
rm -rf ${deck}.${compiler}
clean=YES
if [ ! -d "${deck}.${compiler}" ]; then
  git clone $MODELROOT ${deck}.${compiler}
  clean=NO
fi

cd $REGSCRATCH/${deck}.${compiler}/decks
make rundeck RUN=$deck RUNSRC=$deck >> make.log.${compiler} 2>&1
wait
if [ "\$clean" == "YES" ]; then
  make clean >> make.log.${compiler} 2>&1
  wait
fi
EOF

  if [ "$compiler" == "intel" ]; then

    cat << EOF >> $jobScript
make -j gcm RUN=$deck EXTRA_FFLAGS="-O0 -g -traceback" MPI=$mpi >> make.log.${compiler} 2>&1
wait
EOF

  else

    cat << EOF >> $jobScript
make -j gcm RUN=$deck EXTRA_FFLAGS="-O0 -g -fbacktrace" MPI=$mpi  >> make.log.${compiler} 2>&1
wait
EOF

  fi

  cat << EOF >> $jobScript
make tests RUN=$deck MPI=$mpi >> $testLog 2>&1
wait
EOF
  chmod +x $jobScript

# SUBMIT JOB SCRIPT

  echo "RESULTS [$compiler MPI=$mpi]:" >> $toEmail
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

  local anyError=`grep ' Error' $testLog` 
  if [ "$anyError" != "" ]; then
    local errMsg="Error detected during unit tests" >> $toEmail
    echo "SUMMARY:"  >> $toEmail
    echo "..."  >> $toEmail
    tail -20 $testLog >> $toEmail
    echo ""  >> $toEmail
  fi 
   
  local testFail=`grep 'Failure in' $testLog` 
  if [ "$testFail" != "" ]; then
    nt=`cat $testLog | wc | awk '{print $1}'`
    ni=`cat $testLog | grep -in 'mpirun -np' | awk -F: '{print $1}'`
    n=$((nt-ni)); nh=$((n-1)) 
    # don't want to see make error message, so ...
    tail -$n $testLog > foo; head -$nh foo >> $toEmail; rm foo
  else
    msg=$(tail -2 $testLog | grep " run")
    if [ "$msg" == "" ]; then
      cp $testLog ${testLog}.${pfunitSuffix}.FAILED
      echo "Tests failed to run. Please check ${testLog}.${pfunitSuffix}.FAILED"   >> $toEmail
    else
      tail -2 $testLog | grep " run" >> $toEmail
      rm -f $testLog
    fi
    echo ""  >> $toEmail
  fi 
  wait

# DELETE FILES
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
mpiMode=(YES NO)
for mpi in "${mpiMode[@]}"; do  
  for compiler in "${compilers[@]}"; do 
    job=modelE.${compiler}.j
    log=${ROOT}${compiler}".log"
    submitJob "$compiler" "$job" "$log" "$mpi"
  done 
done 

exit 0
