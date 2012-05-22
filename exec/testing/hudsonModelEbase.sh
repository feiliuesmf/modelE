export WORKSPACE=${WORKSPACE}

# Needed by runModelEBase.sh
export MODELEBASE=/discover/nobackup/ccruz/devel/modelE.clones/base

# the following PBS job will run runModelEBase.sh 
jobID=`qsub /discover/nobackup/ccruz/devel/scripts/runModelEBase.j`

# watch the job
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
if [ -e ${WORKSPACE}/runModelEBase.success ]; then exit 0; else exit 1; fi
