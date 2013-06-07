#!/bin/bash
#PBS -l select=1:ncpus=12
#PBS -l walltime=2:00:00
#PBS -W group_list=s1001
#PBS -N diffrepo
#PBS -j oe
#PBS -V

# This job script is executed by mainRegTests.sh after the regression tests
# have been executed. Its main purpose is to compare the restart results
# using diffreport and generate a DiffReport file emailed to giss-modelE-regression 

# Exit/check codes
FILE_ERR=69
EXIT_ERR=1
OK=0
checkMPI=0
checkSERIAL=0

# -------------------------------------------------------------------
updReport()
# -------------------------------------------------------------------
{
   local message=$1
   echo $message
   line=`echo -e $message`
   report=( "${report[@]}" "$line" )
}

# -------------------------------------------------------------------
checkStatus() 
# -------------------------------------------------------------------
{
   local rc=$1
   local name1=$2
   local name2=$3
   local file1=${name1##*/}
   local file2=${name2##*/}
   if [ $rc -ne $OK ]; then
      if [[ "$name2" =~ baseline ]]; then
         updReport "  ->WARNING: $file1 BASELINE CHANGED"
         deckResults[1]="NO"
      elif [[ "$name1" =~ MPI && "$name2" =~ SERIAL ]]; then
         updReport "  ->ERROR comparing $file1 MPI vs SERIAL"
         deckResults[3]="NO"
      else
         updReport "  ->ERROR in $file1 restart reproducibility"
         deckResults[2]="NO"
      fi
      export willPrintAdditional=YES
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
      updReport "  ->ERROR: $file does NOT exist"
      export willPrintAdditional=YES
      # Model failed during compilation or at runtime
      deckResults[0]="NO"
      deckResults[1]="---"
      deckResults[2]="---"
      deckResults[3]="---"
      return $FILE_ERR
   else
      # Model ran
      deckResults[0]="OK"
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
           if [ $diffSize -ne 0 ]; then
             # Update baseline directory"
             cp -f $name1 $name2
           fi
         fi
         checkStatus $diffSize "$name1" "$name2"
         return $?
      else
         return $EXIT_ERR
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
   declare -a deckResults

   for deck in "${deckArray[@]}"; do
      # defaults
      compileErr=OK
      baseNotChanged=YES
      isRstReprod=YES
      isNPEReprod=YES
      deckResults=($compileErr $baseNotChanged $isRstReprod $isNPEReprod)
      export deckResults
      report=( "${report[@]}" "$deck [$comp] :" )
      echo "  --- DECK = $deck ---"
      # Don't do serial comparisons of C90 and AR5 rundecks
      if [[ "$deck" =~ C90 ]] || [[ "$deck" =~ AR5_CAD ]] || [[ "$deck" =~ tomas ]] || [[ "$deck" =~ amp ]]; then
         echo "  ->Skip SERIAL comparison"
      else
        if [ $checkSERIAL -gt 0 ]; then
# compare SERIAL restart reproducibility
          if [[ ! "$deck" =~ SCM ]]; then
            echo "  ->compare SERIAL restart reproducibility..."
            doDiff $deck.SERIAL.$comp.1dy $deck.SERIAL.$comp.restart $deck $comp
          else
            echo "  ->Skip restart reproducibility..."
          fi
# compare SERIAL baseline (previous day) restart reproducibility
          echo "  ->compare SERIAL baseline reproducibility.."
          doDiff $deck.SERIAL.$comp.1hr $baseline/$deck.SERIAL.$comp.1hr $deck $comp
          doDiff $deck.SERIAL.$comp.1dy $baseline/$deck.SERIAL.$comp.1dy $deck $comp
        fi
      fi
      if [[ "$comp" =~ nag ]] || [[ "$deck" =~ SCM ]]; then
        echo "  ->Skip MPI comparisons when using NAG compiler or SCM rundeck"
      else
# compare MPI restart reproducibility - 3rd argument ($3) is NPE configuration
        if [ ! -z $3 ]; then
          declare -a npeArray=("${!3}")
          echo "  ->compare MPI restart reproducibility..."
          for npe in "${npeArray[@]}"; do
            if [ $checkMPI -gt 0 ]; then
              doDiff $deck.MPI.$comp.1dy.np=$npe $deck.MPI.$comp.restart.np=$npe $deck $comp
            fi
          done
# compare MPI baseline (previous day) restart reproducibility
          echo "  ->compare MPI baseline reproducibility..."
          for npe in "${npeArray[@]}"; do
            if [ $checkMPI -gt 0 ]; then
              doDiff $deck.MPI.$comp.1hr.np=$npe $baseline/$deck.MPI.$comp.1hr.np=$npe $deck $comp
              doDiff $deck.MPI.$comp.1dy.np=$npe $baseline/$deck.MPI.$comp.1dy.np=$npe $deck $comp
            fi
# compare MPI vs SERIAL reproducibility
            echo "  ->compare MPI vs SERIAL reproducibility..."
            if [ $checkMPI -gt 0 ]; then
	      if [[ "$deck" =~ C90 ]] || [[ "$deck" =~ AR5_CAD ]] || [[ "$deck" =~ tomas ]] || [[ "$deck" =~ amp ]] || [[ "$comp" =~ nag ]] || [[ "$deck" =~ SCM ]]; then
                echo "  ->SKIP compare MPI vs SERIAL reproducibility.."
              else
	        doDiff $deck.MPI.$comp.1hr.np=$npe $deck.SERIAL.$comp.1hr $deck $comp
	        doDiff $deck.MPI.$comp.1dy.np=$npe $deck.SERIAL.$comp.1dy $deck $comp
              fi
            fi
          done
        fi
      fi # skip MPI comparisons
      CAS="$deck $comp ${deckResults[@]}"
      deckReport=( "${deckReport[@]}" "$CAS" )
   done
}

# -------------------------------------------------------------------
checkENVS()
# -------------------------------------------------------------------
{
   if [ -z $CONFIG ]; then
      echo " *** ERROR ***"
      echo "ENV variable CONFIG is not defined."
      exit 1;
   else
      echo "CONFIG ENV: $CONFIG"
   fi

   if [ -z $MOCKMODELE ]; then
      diffExec=diffreport.x
      command -v $diffExec &>/dev/null || 
      { 
         echo " ~~~ $diffExec does not exist. Will use GISS installation." >&2
         diffExec=/discover/nobackup/projects/giss/exec/diffreport
      }
   else
      diffExec=/usr/bin/cmp
   fi
}

# -------------------------------------------------------------------
readCFG()
# -------------------------------------------------------------------
{
# save IFS = internal field separator (default is space/tab/newline)
   OIFS=$IFS
# "=" separates fields
   IFS="="

# cfg is the file created by the regresssionTests.pl script
   cfg="$TESTD/."$CONFIG".cfg"
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
# restore IFS
   IFS=$OIFS

   echo "Total DECKs: ${id}"
}

# -------------------------------------------------------------------
separateDecks()
# -------------------------------------------------------------------
{
# NPES used varies by rundeck depending on resolution
# Unfortunately the config file does not contain this information, which
# is part of the file names. Therefore, it must be specified here and must
# consistent with the values in regressionTests.pl.
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
   elif [[ "$LEVEL" == "XLDECK" ]]; then
      LowResNpes=( )
      HiResNpes=( 44 )
      CSNpes=( )
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
      if [[ "$deck" =~ EM20 || "$deck" =~ E1oM20  || "$deck" =~ C12 ]]; then
        LowResDecks[$ia]="$deck"
        ia=$(($ia+1))
      elif [[ "$deck" =~ obio || "$deck" =~ cadF40 || "$deck" =~ amp || "$deck" =~ tomas ]]; then
        HiResDecks[$ib]="$deck"
        ib=$(($ib+1))
      elif [[ "$deck" =~ C90  ]]; then
        CSDecks[$ic]="$deck"
        ic=$(($ic+1))
      elif [[ "$deck" =~ AR5_CAD ]]; then
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
}

# -------------------------------------------------------------------
createEmailReport()
# -------------------------------------------------------------------
{
# Create report for email

   rm -f $TESTD/${CONFIG}.diff
   echo "ModelE test results, branch=$branch" 
   echo "--------------------------------------------------------------------------"

   len=${#deckReport[*]}
   i=0
   while [ $i -lt $len ]; do
      echo "${deckReport[$i]}" 
      echo "${deckReport[$i]}" >> $TESTD/.diffrep
      let i++
   done

echo "ModelE test results, branch=$branch" > $TESTD/.foo
echo "------------------------------------------------------------" >> $TESTD/.foo
echo "                                  -REPRODUCIBILITY-" >> $TESTD/.foo
echo "        RUNDECK""   COMPILER ""  RUN ""  BAS ""  RST ""  NPE " >> $TESTD/.foo
echo "------------------------------------------------------------" >> $TESTD/.foo
awk '{ printf "%15s %10s %5s %5s %5s %5s\n", $1, $2, $3, $4, $5, $6 }' $TESTD/.diffrep  >> $TESTD/.foo

   mv $TESTD/.foo $TESTD/${CONFIG}.diff
   if [ "$willPrintAdditional" == "YES" ]; then
      for ((i=0; i < ${#report[@]}; i++)); do 
         echo "${report[${i}]}" 
         #echo "${report[${i}]}" >> $TESTD/${CONFIG}.diff
      done
   fi
}

# -------------------------------------------------------------------
copy2Workspace()
# -------------------------------------------------------------------
{
# HUDSON monitors a workspace where it can access copied files
   if [ -z $WORKSPACE ]; then
      echo " *** WARNING ***"
      echo "WORKSPACE is not defined. HUDSON will report a failure."
   else
      echo "Imported WORKSPACE = "$WORKSPACE
      rm -f $WORKSPACE/.success
      writeOK=1
   fi
# copy the generated email reports
   cp -f $TESTD/${CONFIG}.diff $WORKSPACE
   cp -f $TESTD/${CONFIG}.unit $WORKSPACE
}

# -------------------------------------------------------------------
archive()
# -------------------------------------------------------------------
{
# Archive full difference reports
   cp -f $TESTD/${CONFIG}.diff $MODELEBASELINE/reports/${CONFIG}.diff.`date +%F`

# Check for errors in report
   cat $TESTD/${CONFIG}.diff | grep -i " NO " > /dev/null
   rc1=$?
   cat $TESTD/${CONFIG}.diff | grep -i "ERR" > /dev/null
   rc2=$?
# if we found errors, (0 return code!) then we exit
   if [ $rc1 -eq 0 ] || [ $rc2 -eq 0 ]; then
      echo "Regression tests ERRORS: Will NOT create modelE snapshot"
      exit $EXIT_ERR
   else 
# Create modelE snapshot iff no ERRORs in ${CONFIG}.diff (WARNINGs are OK)
      if [ "$writeOK" -eq 1 ]; then touch $WORKSPACE/.success; fi
      echo "Will create modelE snapshot"
      # Create modelE snapshot
      if [ -z $MOCKMODELE ]; then
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
   fi
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

export willPrintAdditional=NO

TESTD=$MODELROOT/exec/testing

cd $TESTD

rm -f .diffrep

checkENVS

# read config file created by regressionTests.pl
readCFG

# organize decks into sets identified by size and NPES used
separateDecks

# Loop over cases and diff the results
for comp in "${COMPILERS[@]}"; do

  echo "--- COMPILER = $comp ---"

  cd $REGRESULTS/$comp

  checkMPI=`ls -1 *MPI* | wc -c`
  checkSERIAL=`ls -1 *SERIAL* | wc -c`

  report=( "${report[@]}" "" )
  report=( "${report[@]}" "ADDITIONAL DETAILS:")
  report=( "${report[@]}" "===================")
  deckDiff $comp LowResDecks[@] LowResNpes[@]
  deckDiff $comp HiResDecks[@] HiResNpes[@]
  deckDiff $comp CSDecks[@] CSNpes[@]
  deckDiff $comp AR5Decks[@] HiResNpes[@]
  deckDiff $comp SCMdecks[@]

done

createEmailReport

copy2Workspace

archive

exit $OK
