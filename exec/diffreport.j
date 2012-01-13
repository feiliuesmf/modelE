#!/bin/bash
#PBS -l select=1
#PBS -l walltime=2:00:00
#PBS -W group_list=a940a
#PBS -N diffrep
#PBS -j oe
#PBS -V

# This job script is executed by mainRegTests.sh after the regression tests
# have been executed. Its main purppose is to compare the restart results
# using diffreport and generates a DiffReport file emailed to giss-modelE-regression 
# Note:
# This script is a temporary measure. It has a couple of "dependencies" and various
# hardwired directory paths. Also, it would be highly desirable to make the
# restart comparisons within the perl regression test scripts but somehow that
# task takes several hours to complete. Running within this job takes ~30mins.

cd $PBS_O_WORKDIR
rm -f DiffReport diffrep.o*

# Dependency: I use /home/modele/exec/diffreport.x 
diffDiffReport=/home/modele/exec/diffreport.x
declare -a report
declare -a DECKS
declare -a COMPILERS
declare -a LowResDecks
declare -a HiResDecks
declare -a SCMdecks

OIFS=$IFS
IFS="="

# Read configuration file...always stored in master/ directory
cfg=/home/modele/master/exec/.regTest.cfg
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
      elif [[ "$deck" =~ E4F40 || "$deck" =~ E4Tcad || "$deck" =~ E4arobio ]]; then
        HiResDecks[$ib]="$deck"
        ib=$(($ib+1))
      else
        SCMdecks[$ic]="$deck"
        ic=$(($ic+1))
      fi
done

echo "LowResDecks: ${LowResDecks[@]}"
echo "HiResDecks: ${HiResDecks[@]}"
echo "SCMdecks: ${SCMdecks[@]}"
echo "COMPILERS: ${COMPILERS[*]}"
echo "LEVEL is: $LEVEL"
echo "BRANCH is: $BRANCH"

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
   fileExists "$name1"
   export isReprod=YES
   if [ $? -eq 0 ]; then
      fileExists "$name2"
      if [ $? -eq 0 ]; then
         $diffDiffReport $name1 $name2 > foo
         wait
         foosize=`cat foo | wc -c`; rm -f foo
         # save file to BASELINE directory
         if [[ "$name2" =~ baseline ]]; then
         #   if [ $foosize -eq 0 ]; then
               echo " Updating $deck baseline"
               cp -f $name1 $name2
         #   fi
         fi
         checkStatus $foosize "$name1" "$name2"
         return $?
      fi
   fi
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


for comp in "${COMPILERS[@]}"; do

  echo "--- COMPILER = $comp ---"

  cd /discover/nobackup/modele/regression_results/$comp
  baseline=/discover/nobackup/modele/modelE_baseline/$comp

  for deck in "${LowResDecks[@]}"; do

  echo "--- DECK = $deck ---"

# compare SERIAL restart reproducibility
    doDiff $deck.SERIAL.$comp.1dy $deck.SERIAL.$comp.restart $deck
# compare MPI restart reproducibility
    for npe in "${LowResNpes[@]}"; do
      doDiff $deck.MPI.$comp.1dy.np=$npe $deck.MPI.$comp.restart.np=$npe $deck
    done
# compare SERIAL baseline (previous day) restart reproducibility
    doDiff $deck.SERIAL.$comp.1hr $baseline/$deck.SERIAL.$comp.1hr $deck
    doDiff $deck.SERIAL.$comp.1dy $baseline/$deck.SERIAL.$comp.1dy $deck
# compare MPI baseline (previous day) restart reproducibility
    for npe in "${LowResNpes[@]}"; do
      doDiff $deck.MPI.$comp.1hr.np=$npe $baseline/$deck.MPI.$comp.1hr.np=$npe $deck
      doDiff $deck.MPI.$comp.1dy.np=$npe $baseline/$deck.MPI.$comp.1dy.np=$npe $deck
    done
    if [ "$isReprod" == "YES" ]; then
       upd "$deck" "$comp" "STRONGLY" 
    else
       upd "$deck" "$comp" "NOT" 
    fi

  done

  for deck in "${HiResDecks[@]}"; do

  echo "--- DECK = $deck ---"

# compare SERIAL restart reproducibility
    doDiff $deck.SERIAL.$comp.1dy $deck.SERIAL.$comp.restart $deck
# compare MPI restart reproducibility
    for npe in "${HiResNpes[@]}"; do
      doDiff $deck.MPI.$comp.1dy.np=$npe $deck.MPI.$comp.restart.np=$npe $deck
    done
# compare SERIAL baseline (previous day) restart reproducibility
    doDiff $deck.SERIAL.$comp.1hr $baseline/$deck.SERIAL.$comp.1hr $deck
    doDiff $deck.SERIAL.$comp.1dy $baseline/$deck.SERIAL.$comp.1dy $deck
# compare MPI baseline (previous day) restart reproducibility
    for npe in "${HiResNpes[@]}"; do
      doDiff $deck.MPI.$comp.1hr.np=$npe $baseline/$deck.MPI.$comp.1hr.np=$npe $deck
      doDiff $deck.MPI.$comp.1dy.np=$npe $baseline/$deck.MPI.$comp.1dy.np=$npe $deck
    done
    if [ "$isReprod" == "YES" ]; then
       upd "$deck" "$comp" "STRONGLY" 
    else
       upd "$deck" "$comp" "NOT" 
    fi

  done

  for deck in "${SCMdecks[@]}"; do

  echo "--- DECK = $deck ---"

# compare SERIAL restart reproducibility
    doDiff $deck.SERIAL.$comp.1dy $deck.SERIAL.$comp.restart $deck
# compare SERIAL baseline (previous day) restart reproducibility
#    doDiff $deck.SERIAL.$comp.1hr $baseline/$deck.SERIAL.$comp.1hr $deck
#    doDiff $deck.SERIAL.$comp.1dy $baseline/$deck.SERIAL.$comp.1dy $deck
    if [ "$isReprod" == "YES" ]; then
       upd "$deck" "$comp" "STRONGLY" 
    else
       upd "$deck" "$comp" "NOT" 
    fi

  done

done

touch $HOME/master/exec/DiffReport
echo "Results:"
for ((i=0; i < ${#report[@]}; i++)); do 
   echo "${report[${i}]}"
   echo "${report[${i}]}" >> $HOME/master/exec/DiffReport
done

cat $HOME/master/exec/DiffReport | grep ERROR > /dev/null 2>&1
RC=$?
# Create modeleE snapshot iff no ERRORs in DiffReport
if [ $RC -eq 0 ]; then
   echo "Regression tests ERROR: Will NOT create modelE snapshot"
else 
   echo "Will create modelE snapshot"
   # Create modeleE snapshot
   if [ -d "$NOBACKUP/regression_scratch/modelE" ]; then
      cd $NOBACKUP/regression_scratch/$BRANCH
      DST=/discover/nobackup/modele/modelE_baseline/snapshots/
      NAME=modelE.`date +%F`.zip
      git archive -o $DST/$NAME $BRANCH
   else 
      ls $NOBACKUP/regression_scratch/$BRANCH
      echo "Could not create modelE snapshot"
   fi
fi

exit 0
