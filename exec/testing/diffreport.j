#!/usr/local/bin/bash
#PBS -l select=1:ncpus=12:mpiprocs=12
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

# Exit codes
FILE_ERR=69
EXIT_ERR=1
OK=0

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
   if [ $rc -ne $OK ]; then
      if [ $rc -gt 1 ]; then
         updDiffReport " --- ERROR: One or more missing files"
      else
         if [[ "$name2" =~ baseline ]]; then
            updDiffReport " --- WARNING: inconsistent results with BASELINE $base1"
         else
            updDiffReport " --- WARNING: inconsistent results between $base1 and $base2"
         fi
      fi
      export isReprod=NO
      return $E_EXIT_ERR
   fi
   return $OK
}

# -------------------------------------------------------------------
fileExists()
# -------------------------------------------------------------------
{
   local file=$1
   #local file=${1##*/}
   if [ ! -e "$file" ]; then
      updDiffReport " --- ERROR: $file does NOT exist"
      export isReprod=NO
      return $FILE_ERR
   fi
   return $OK
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
   return_val=$?
   if [ "$return_val" -eq $OK ]; then
      fileExists "$name2"
      return_val=$?
      if [ "$return_val" -eq $OK ]; then
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
      return $EXIT_ERR
   fi
}

# -------------------------------------------------------------------
deckDiff()
# -------------------------------------------------------------------
{
  local comp="$1"
  declare -a deckArray=("${!2}")
  # if deckArray is empty then there is nothing to do:
  if [ ${#deckArray[@]} -eq 0 ]; then return; fi
  local baseline=$NOBACKUP/modelE_baseline/$comp

  for deck in "${deckArray[@]}"; do

    echo "  --- DECK = $deck ---"
    # Don't do serial comparisons of C90 and AR5 rundecks
    if [[ "$deck" =~ C90 ]] || [[ "$deck" =~ AR5 ]]; then
      echo "    Skip SERIAL comparison"
    else
# compare SERIAL restart reproducibility
      doDiff $deck.SERIAL.$comp.1dy $deck.SERIAL.$comp.restart $deck $comp
# compare SERIAL baseline (previous day) restart reproducibility
      doDiff $deck.SERIAL.$comp.1hr $baseline/$deck.SERIAL.$comp.1hr $deck $comp
      doDiff $deck.SERIAL.$comp.1dy $baseline/$deck.SERIAL.$comp.1dy $deck $comp
    fi
# compare MPI restart reproducibility - 3rd argument ($3) is NPE configuration
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
    if [ "$isReprod" == YES ]; then
      upd "$deck" "$comp" "REPRODUCIBLE" 
    else
      upd "$deck" "$comp" "***NOT REPRODUCIBLE***" 
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
  updDiffReport " $deck_ [ $comp_ ] is $which_" 
}

# -------------------------------------------------------------------
# MAIN
# -------------------------------------------------------------------

cd $PBS_O_WORKDIR
rm -f diffrep.o*

if [ -z $MOCKMODELE ]; then
  diffDiffReport=$NOBACKUP/devel/master/exec/testing/diffreport.x
  if [ ! -e $diffDiffReport ]; then
     echo " *** WARNING ***"
     echo "$diffDiffReport does not exist"
     echo "Will use unix cmp but Diffreport may be inaccurate"
     diffDiffReport=/usr/bin/cmp
  fi
else
  diffDiffReport=/usr/bin/cmp
fi

declare -a report
declare -a DECKS
declare -a COMPILERS
declare -a LowResDecks
declare -a HiResDecks
declare -a CSDecks
declare -a AR5Decks
declare -a SCMdecks

OIFS=$IFS
IFS="="

# Read configuration file...always stored in master/ directory
cfg=$NOBACKUP/devel/master/exec/testing/.regTest.cfg
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

# Also NPES used varies by rundeck depending on resolution noted above
if [[ "$LEVEL" == "GENTLE" ]]; then
   LowResNpes=( 4 )
   HiResNpes=( 8 )
   CSNpes=( 6 )
elif [[ "$LEVEL" == "AGGRESSIVE" ]]; then
   LowResNpes=( 1 4 23 )
   HiResNpes=( 1 45 )
   CSNpes=( 48 )
elif [[ "$LEVEL" == "INSANE" ]]; then
   LowResNpes=( 1 4 23 44 )
   HiResNpes=( 1 8 45 88 )
   CSNpes=( 84 )
else
   LowResNpes=( 4 )
   HiResNpes=( 8 )
   CSNpes=( 6 )
fi

# Separate rundecks to test into 4x4.5, 2x2.5 and SCM
ia=0
ib=0
ic=0
ir=0
is=0
for deck in ${DECKS[@]}
do
      if [[ "$deck" =~ EM20 || "$deck" =~ E1oM20 ]]; then
        LowResDecks[$ia]="$deck"
        ia=$(($ia+1))
      elif [[ "$deck" =~ obio || "$deck" =~ cadF40  ]]; then
        HiResDecks[$ib]="$deck"
        ib=$(($ib+1))
      elif [[ "$deck" =~ C90  ]]; then
        CSDecks[$ic]="$deck"
        ic=$(($ic+1))
      elif [[ "$deck" =~ AR5  ]]; then
        AR5Decks[$ir]="$deck"
        ir=$(($ir+1))
      else
        SCMdecks[$is]="$deck"
        is=$(($is+1))
      fi
done

echo "LowResDecks: ${LowResDecks[@]}"
echo "HiResDecks: ${HiResDecks[@]}"
echo "SCMdecks: ${SCMdecks[@]}"
echo "CSDecks: ${CSDecks[@]}"
echo "AR5Decks: ${AR5Decks[@]}"
echo "COMPILERS: ${COMPILERS[*]}"
echo "LEVEL is: $LEVEL"
echo "BRANCH is: $BRANCH"
echo "LowResNpes is: ${LowResNpes[@]}"
echo "HiResNpes is: ${HiResNpes[@]}"
echo "CSNpes is: ${CSNpes[@]}"

for comp in "${COMPILERS[@]}"; do

  echo "--- COMPILER = $comp ---"

  cd $NOBACKUP/regression_results/$comp

  deckDiff $comp LowResDecks[@] LowResNpes[@]
  deckDiff $comp HiResDecks[@] HiResNpes[@]
  deckDiff $comp CSDecks[@] CSNpes[@]
  deckDiff $comp AR5Decks[@] HiResNpes[@]
  deckDiff $comp SCMdecks[@]

done

rm -f $NOBACKUP/devel/master/exec/testing/DiffReport
echo "Results:"
for ((i=0; i < ${#report[@]}; i++)); do 
   echo "${report[${i}]}"
   echo "${report[${i}]}" >> $NOBACKUP/devel/master/exec/testing/DiffReport
done

cat $NOBACKUP/devel/master/exec/testing/DiffReport | grep ERROR > /dev/null 2>&1
rc=$?
# Create modelE snapshot iff no ERRORs in DiffReport (WARNINGs are OK)
if [ $rc -eq $OK ]; then
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

exit $OK
