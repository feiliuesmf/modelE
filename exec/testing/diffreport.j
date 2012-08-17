#!/bin/bash
#PBS -l select=1:ncpus=12:mpiprocs=12
#PBS -l walltime=2:00:00
#PBS -W group_list=s1001
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
updReport()
# -------------------------------------------------------------------
{
   local message=$1
   line=`echo -e $message`
   report=( "${report[@]}" "$line" )
}

# -------------------------------------------------------------------
updDeckReport()
# -------------------------------------------------------------------
{
   local message=$1
   line=$(echo $message | awk '{ printf "%-15s | %-10s | %-30s\n", $1,$2,$3}')
   deckReport=( "${deckReport[@]}" "$line" )
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
      if [[ "$name2" =~ baseline ]]; then
         echo " --- WARNING: $base1 BASELINE has CHANGED"
         updReport " --- WARNING: $base1 BASELINE has CHANGED"
	 export isChanged=YES
      else
         echo " --- ERROR: $base1 and $base2 DIFFER"
         updReport " --- ERROR: $base1 and $base2 DIFFER"
	 export isReprod=NO
      fi
      export printReport=YES
      return $E_EXIT_ERR
   fi
   return $OK
}

# -------------------------------------------------------------------
fileExists()
# -------------------------------------------------------------------
{
   local file=$1
   if [ ! -e "$file" ]; then
      updReport " --- ERROR: $file does NOT exist"
      export compileErr=YES
      export printReport=YES
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
   fileExists "$name1"
   return_val=$?
   if [ "$return_val" -eq $OK ]; then
      fileExists "$name2"
      return_val=$?
      if [ "$return_val" -eq $OK ]; then
         $diffExec $name1 $name2 > fileDiff
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
             # save old BASELINE 
             cp -f $name2 $name2.save
             # replace BASELINE 
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
  local baseline=$MODELEBASELINE/$comp

  for deck in "${deckArray[@]}"; do
    export compileErr=NO
    export isChanged=NO
    export isReprod=YES
    report=( "${report[@]}" "$deck [$comp] :" )
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
    if [ "$compileErr" == YES ]; then
      updDeckReport "$deck $comp ***COMPILE_RUNTIME_ERROR***"
    else
      if [ "$isReprod" == YES ]; then
         updDeckReport "$deck $comp IS_reproducible"
         if [ "$isChanged" == YES ]; then
            updDeckReport "$deck $comp CHANGED_baseline"
         fi
      else
         updDeckReport "$deck $comp NOT_reproducible"
      fi
    fi
  done

}

# -------------------------------------------------------------------
# MAIN
# -------------------------------------------------------------------
umask 002

declare -a deckReport
declare -a report
declare -a DECKS
declare -a COMPILERS
declare -a LowResDecks
declare -a HiResDecks
declare -a CSDecks
declare -a AR5Decks
declare -a SCMdecks

export printReport=NO

cd $MODELROOT/exec/testing

if [ -z $CONFIG ]; then
   echo " *** ERROR ***"
   echo "ENV variable CONFIG is not defined."
   exit 1;
else
   echo "CONFIG ENV: $CONFIG"
fi

if [ -z $MOCKMODELE ]; then
  compiler=$(cat $MODELERC | grep COMPILER= | awk -F"=" '{print $2}')
  if [ "$compiler" == "intel" ]; then
     FC=ifort
  else
     FC=gfortran
  fi
  netcdf=$(cat $MODELERC | grep "/netcdf/" | awk -F"=" '{print $2}')
  INC=$netcdf/include
  LIB=$netcdf/lib

  misnc=$MODELROOT/model/mk_diags/miscnc.f
  diffr=$MODELROOT/model/mk_diags/diffreport.f
  diffExec=$MODELROOT/exec/testing/diffreport.x
  
  echo $FC -O3 -I$INC $misnc $diffr -o $diffExec -L$LIB -lnetcdf
  $FC -O3 -I$INC $misnc $diffr -o $diffExec -L$LIB -lnetcdf
  if [ ! -e $diffExec ]; then
     echo " *** WARNING ***"
     echo "$diffExec does not exist. Will use GISS installation"
     diffExec=/discover/nobackup/projects/giss/exec/diffreport
  fi
else
  diffExec=/usr/bin/cmp
fi

OIFS=$IFS
IFS="="

cfg="$MODELROOT/exec/testing/."$CONFIG".cfg"
if [ ! -e $cfg ]; then
   echo " *** ERROR ***"
   echo "$cfg does not exist."
   exit 1
else
   echo "Config file: $cfg"
fi

id=0
ic=0
# Read configuration file $cfg
while read line ; do
   set -- $line
   arr=($*)
   if [[ "${arr[0]}" == "BRANCH" ]]; then
      branch=${arr[1]}
      BRANCH=${arr[1]}
   fi
   if [[ "${arr[0]}" == "LEVEL" ]]; then
      level=${arr[1]}
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

  cd $REGRESULTS/$comp

  report=( "${report[@]}" "" )
  report=( "${report[@]}" "ADDITIONAL DETAILS:")
  report=( "${report[@]}" "===================")
  deckDiff $comp LowResDecks[@] LowResNpes[@]
  deckDiff $comp HiResDecks[@] HiResNpes[@]
  deckDiff $comp CSDecks[@] CSNpes[@]
  deckDiff $comp AR5Decks[@] HiResNpes[@]
  deckDiff $comp SCMdecks[@]

done

rm -f $MODELROOT/exec/testing/${CONFIG}.diff
echo "ModelE test results, branch=$branch" 
echo "--------------------------------------------------------------------------"
echo "ModelE test results, branch=$branch" >> $MODELROOT/exec/testing/${CONFIG}.diff
echo "--------------------------------------------------------------------------" >> $MODELROOT/exec/testing/${CONFIG}.diff

len=${#deckReport[*]}
i=0
while [ $i -lt $len ]; do
   echo "${deckReport[$i]}" 
   echo "${deckReport[$i]}" >> $MODELROOT/exec/testing/${CONFIG}.diff
   let i++
done

if [ "$printReport" == "YES" ]; then
   for ((i=0; i < ${#report[@]}; i++)); do 
      echo "${report[${i}]}" 
      echo "${report[${i}]}" >> $MODELROOT/exec/testing/${CONFIG}.diff
   done
fi

if [ -z $WORKSPACE ]; then
   echo " *** WARNING ***"
   echo "WORKSPACE is not defined. HUDSON will report a failure."
else
   echo "Imported WORKSPACE = "$WORKSPACE
   rm -f $WORKSPACE/.success
   writeOK=1
fi

#chmod g+rw $MODELROOT/exec/testing/${CONFIG}.diff
cp $MODELROOT/exec/testing/${CONFIG}.diff $WORKSPACE
# Archive full difference reports
cp $MODELROOT/exec/testing/${CONFIG}.diff $MODELEBASELINE/reports/${CONFIG}.diff.`date +%F`

# Check for errors in report
cat $MODELROOT/exec/testing/${CONFIG}.diff | grep -i "NOT_REPRODUCIBLE" > /dev/null
rc1=$?
cat $MODELROOT/exec/testing/${CONFIG}.diff | grep -i "ERROR" > /dev/null
rc2=$?
# if we found errors, (0 return code!) then we exit
if [ $rc1 -eq 0 ] || [ $rc2 -eq 0 ]; then
   echo "Regression tests ERRORS: Will NOT create modelE snapshot"
   exit $EXIT_ERR
else 
# Create modelE snapshot iff no ERRORs in ${CONFIG}.diff (WARNINGs are OK)
   if [ $writeOK -eq 1 ]; then touch $WORKSPACE/.success; fi
   echo "Will create modelE snapshot"
   # Create modelE snapshot
   if [ -d "$REGSCRATCH/$BRANCH" ]; then
      cd $REGSCRATCH/$BRANCH
      DST=$MODELEBASELINE/snapshots/
      NAME=modelE.`date +%F`.zip
      git archive -o $DST/$NAME $BRANCH
   else 
      ls $REGSCRATCH/$BRANCH
      echo "Could not create modelE snapshot"
   fi
fi

exit $OK
