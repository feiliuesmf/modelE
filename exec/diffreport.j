#!/bin/bash
#PBS -l select=1
#PBS -l walltime=2:00:00
#PBS -W group_list=a940a
#PBS -N diffrep
#PBS -j oe
#PBS -V

# This job script is executed by mainRegTests.sh after the regression tests
# have been executed. Its main purpose is to compare the restart results
# using diffreport and generate a DiffReport file emailed to giss-modelE-regression 
# Note:
# This script is a temporary measure. It has a couple of "dependencies" and various
# hardwired directory paths. Also, it would be highly desirable to make the
# restart comparisons within the perl regression test scripts but somehow that
# task takes several hours to complete. Running within this job takes ~30mins.

# -------------------------------------------------------------------
updDiffReport()
# -------------------------------------------------------------------
{
   local message=$1
   line=`echo -e $message"\n"`
   report=( "${report[@]}" "$line" )
}

# -------------------------------------------------------------------
checkStatus() 
# -------------------------------------------------------------------
{
   local rc=$1
   local name1=$2
   local name2=$3
   local base1=${name1##*/}
   local base2=${name2##*/}
   if [ $rc -ne 0 ]; then
      if [[ "$name2" =~ baseline ]]; then
         updDiffReport " --- WARNING: inconsistent results with BASELINE $base1"
      else
         updDiffReport " --- WARNING: inconsistent results between $base1 and $base2"
      fi
      export isReprod=NO
      return 1
   fi
   return 0
}

# -------------------------------------------------------------------
fileExists()
# -------------------------------------------------------------------
{
   local file=$1
   if [ ! -e "$file" ]; then
      updDiffReport " ERROR: $file does NOT exist"
      export isReprod=NO
      return 1
   fi
   return 0
}

# -------------------------------------------------------------------
doDiff()
# -------------------------------------------------------------------
{ 
   local name1=$1
   local name2=$2
   local deck=$3
   local comp=$4
   export isReprod=YES
   fileExists "$name1"
   if [ $? -eq 0 ]; then
      fileExists "$name2"
      if [ $? -eq 0 ]; then
         $diffDiffReport $name1 $name2 > fileDiff
         wait
         diffSize=`cat fileDiff | wc -c`; rm -f fileDiff
         # save file to BASELINE directory
         if [[ "$name2" =~ baseline ]]; then
           if [ $diffSize -eq 0 ]; then
             echo "    $name1 results and baseline are IDENTICAL."
             echo "    -- Will NOT update baseline directory."
           else
             echo "    $name1 results and baseline DIFFER."
             echo "    -- Updating $deck in baseline directory"
             cp -f $name1 $name2
           fi
         fi
         checkStatus $diffSize "$name1" "$name2"
         return $?
      fi
   else
      return 1
   fi
}

# -------------------------------------------------------------------
deckDiff()
# -------------------------------------------------------------------
{
  local comp="$1"
  declare -a deckArray=("${!2}")
  local baseline=$NOBACKUP/modelE_baseline/$comp

  for deck in "${deckArray[@]}"; do

    echo "  --- DECK = $deck ---"
    # Don't do serial comparisions of C() and AR5 rundecks
    if [[ "$deck" =~ C90 ]] || [[ "$deck" =~ AR5 ]]; then
      echo "    Skip SERIAL comparison"
    else
# compare SERIAL restart reproducibility
      doDiff $deck.SERIAL.$comp.1dy $deck.SERIAL.$comp.restart $deck $comp
# compare SERIAL baseline (previous day) restart reproducibility
      doDiff $deck.SERIAL.$comp.1hr $baseline/$deck.SERIAL.$comp.1hr $deck $comp
      doDiff $deck.SERIAL.$comp.1dy $baseline/$deck.SERIAL.$comp.1dy $deck $comp
    fi
# compare MPI restart reproducibility
    if [ ! -z $3 ]; then
      declare -a npeArray=("${!3}")
      for npe in "${npeArray[@]}"; do
        doDiff $deck.MPI.$comp.1dy.np=$npe $deck.MPI.$comp.restart.np=$npe $deck $comp
      done
# compare MPI baseline (previous day) restart reproducibility
      for npe in "${npeArray[@]}"; do
        doDiff $deck.MPI.$comp.1hr.np=$npe $baseline/$deck.MPI.$comp.1hr.np=$npe $deck $comp
        doDiff $deck.MPI.$comp.1dy.np=$npe $baseline/$deck.MPI.$comp.1dy.np=$npe $deck $comp
      done
    fi
    if [ $? -eq 0 ]; then
      upd "$deck" "$comp" "STRONGLY" 
    else
      upd "$deck" "$comp" "NOT" 
    fi

  done
}

# -------------------------------------------------------------------
upd()
# -------------------------------------------------------------------
{
  local deck_=$1
  local comp_=$2
  local which_=$3
  updDiffReport " Rundeck $deck_ with compiler $comp_ is $which_ reproducible" 
}

# -------------------------------------------------------------------
# MAIN
# -------------------------------------------------------------------

cd $PBS_O_WORKDIR
rm -f DiffReport diffrep.o*

if [ -z $MOCKMODELE ]; then
  diffDiffReport=$NOBACKUP/devel/master/exec/diffreport.x
else
  diffDiffReport=/usr/bin/cmp
fi

declare -a report
declare -a DECKS
declare -a COMPILERS
declare -a LowResDecks
declare -a HiResDecks
declare -a SCMdecks

OIFS=$IFS
IFS="="

# Read configuration file...always stored in master/ directory
cfg=$NOBACKUP/devel/master/exec/.regTest.cfg
id=0
ic=0
while read line ; do
   set -- $line
   arr=($*)
   if [[ "${arr[0]}" == "BRANCH" ]]; then
      BRANCH=${arr[1]}
   fi
   if [[ "${arr[0]}" == "LEVEL" ]]; then
      LEVEL=${arr[1]}
   fi
   if [[ "${arr[0]}" == "COMPILER" ]]; then
      COMPILERS[$ic]=${arr[1]}
      ic=$(($ic+1))
   fi
   if [[ "${arr[0]}" == "DECK" ]]; then
      DECKS[$id]="${arr[1]}"
      id=$(($id+1))
   fi
done < $cfg
IFS=$OIFS

echo "Total DECKs: ${id}"

# Separate rundecks to test into 4x4.5, 2x2.5 and SCM
ia=0
ib=0
ic=0
for deck in ${DECKS[@]}
do
      if [[ "$deck" =~ EM20 || "$deck" =~ E1oM20 ]]; then
        LowResDecks[$ia]="$deck"
        ia=$(($ia+1))
      elif [[ "$deck" =~ E4 ]]; then
        HiResDecks[$ib]="$deck"
        ib=$(($ib+1))
      else
        SCMdecks[$ic]="$deck"
        ic=$(($ic+1))
      fi
done

# Also NPES used varies by rundeck depending on resolution noted above
if [[ "$LEVEL" == "GENTLE" ]]; then
   LowResNpes=( 4 )
   HiResNpes=( 8 )
elif [[ "$LEVEL" == "AGGRESSIVE" ]]; then
   LowResNpes=( 1 4 23 )
   HiResNpes=( 1 45 )
elif [[ "$LEVEL" == "INSANE" ]]; then
   LowResNpes=( 1 4 23 44 )
   HiResNpes=( 1 8 45 88 )
else
   LowResNpes=( 4 )
   HiResNpes=( 8 )
fi

echo "LowResDecks: ${LowResDecks[@]}"
echo "HiResDecks: ${HiResDecks[@]}"
echo "SCMdecks: ${SCMdecks[@]}"
echo "COMPILERS: ${COMPILERS[*]}"
echo "LEVEL is: $LEVEL"
echo "BRANCH is: $BRANCH"
echo "LowResNpes is: ${LowResNpes[@]}"
echo "HiResNpes is: ${HiResNpes[@]}"

for comp in "${COMPILERS[@]}"; do

  echo "--- COMPILER = $comp ---"

  cd $NOBACKUP/regression_results/$comp

  deckDiff $comp LowResDecks[@] LowResNpes[@]
  deckDiff $comp HiResDecks[@] HiResNpes[@]
  deckDiff $comp SCMdecks[@]

done

touch $NOBACKUP/devel/master/exec/DiffReport
echo "Results:"
for ((i=0; i < ${#report[@]}; i++)); do 
   echo "${report[${i}]}"
   echo "${report[${i}]}" >> $NOBACKUP/devel/master/exec/DiffReport
done

cat $NOBACKUP/devel/master/exec/DiffReport | grep ERROR > /dev/null 2>&1
RC=$?
# Create modelE snapshot iff no ERRORs in DiffReport (WARNINGs are OK)
if [ $RC -eq 0 ]; then
   echo "Regression tests ERROR: Will NOT create modelE snapshot"
else 
   echo "Will create modelE snapshot"
   # Create modelE snapshot
   if [ -d "$NOBACKUP/regression_scratch/$BRANCH" ]; then
      cd $NOBACKUP/regression_scratch/$BRANCH
      DST=$NOBACKUP/modelE_baseline/snapshots/
      NAME=modelE.`date +%F`.zip
      git archive -o $DST/$NAME $BRANCH
   else 
      ls $NOBACKUP/regression_scratch/$BRANCH
      echo "Could not create modelE snapshot"
   fi
fi

exit 0
