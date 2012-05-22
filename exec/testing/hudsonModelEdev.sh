export WORKSPACE=${WORKSPACE}

# ENV variables used by runModelE.sh
export MODELEBASE=/discover/nobackup/ccruz/devel/modelE.clones/base
export MODELEGIT=/discover/nobackup/ccruz/devel/modelE.clones/dev

# Submit the PBS job
jobID=`qsub /discover/nobackup/ccruz/devel/scripts/runModelE.j`

# Watch the job
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
if [ -e ${WORKSPACE}/runModelE.success ]; then exit 0; else exit 1; fi
