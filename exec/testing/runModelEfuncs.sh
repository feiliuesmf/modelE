# -------------------------------------------------------------------
# FUNCTIONS
# -------------------------------------------------------------------

# -------------------------------------------------------------------
usage() 
# -------------------------------------------------------------------
{
less <<EOF
NAME
     $0 - Script to compile/run modelE rundecks

DESCRIPTION

     Use this script to compile and/or run modelE rundecks through command line options.
     The script compiles/runs within a working git repository as specified in the MODELEGIT 
     environment variable. On DISCOVER it is assumed that the script is executed within a PBS
     interactive session. On non-DISCOVER systems a fortran compiler (gfortran or Intel) 
     and MPI must be available.     
     To use this script to create a modelE baseline set the MODELEBASE environment variable
     and pass the appropriate command line flags.
    

SYNOPSIS

     $0 <...arguments...>  [..options...]

REQUIREMENTS

     MODELEGIT      : environment variable defining git local repository path
                      e.g. export MODELEGIT=<path to git repository>

     MODELEBASE      : environment variable defining a baseline repository path
                      e.g. export MODELEBASE=<path to git repository>

ARGUMENTS

     -h | --help     : this help screen

OPTIONS

     -b | --baseline : create modelE baseline result (default: NO)
     -v | --verbose  : verbose (default: YES)
     -f | --fflags   : use debug FFLAGS to compile (default: YES)
     -c | --clean    : run "make clean" (default: NO) 
     -r | --regress  : perform restart regression (default: NO)
     -x | --check    : check run results against baseline repository (default: NO)
                       If yes, set export MODELEBASE=<path to baseline repository>
     -m | --modelErc : specified location of .modelErc (default: $HOME/.modelErc)
                       OR set export MODELERC=<path to file>
     -d | --deck     : rundeck name (default: EM20)
                       OR set export MODELEDECK=<deck name>

     Additionally, on DISCOVER:

     -e | --modEnv   : module environment (default: intel13)
                             intel13  : comp/intel-13.0.0.079 
                                        mpi/impi-3.2.2.006
                             gcc47    : other/comp/gcc-4.7-20120331
					other/mpi/mvapich2-1.8a2/gcc-4.7-20120331
EXAMPLE

     This screen:
        runModelE.sh -h 
     Create baseline:
        runModelE.sh -b YES
     Do not check agains baseline and perform a make clean
        runModelE.sh -x NO -c YES
     Use gcc47 modules, E4TcadF40 rundeck and the specified modelErc file
        runModelE.sh -e gcc47 -d E4TcadF40 --modelErc=$HOME/.modelErc.G47

AUTHOR

     Carlos Cruz (carlos.a.cruz@nasa.gov), NASA Code 610.3

EOF
}

# -------------------------------------------------------------------
defaults()
# -------------------------------------------------------------------
{
# Defaults. Most of these can be overridden via command-line arguments.
   clear

# Debug compilation flags
   fflags="YES"

# .modelErc file path - the default file is automatically created here
   rcFile=$HOME/.modelErc

# Rundeck name
   rundeck="nonProduction_E_AR5_C12"

# module environment in .modelErc - default is to use Intel 13
   moduleSet="intel13"
   # used in .modelErc:
   COMPILER=intel
   MPIDISTR=intel
   MPIDIR=

# Run make clean
   clean="NO"

# Default: Do not do restart regression
   regress="NO"

# Default: Do not create a baseline
   baseline="NO"

# Check rundeck changes against those in BASELINE repository
   check="NO"

# Print out run information
   verbose="YES"

# "diff report" for nc files based on GISS utility. On DISCOVER it is 
# installed and ready to go but user must build it on other systems (and
# make sure it is in the PATH).
   if [[ "$node" =~ borg || "$node" =~ discover || "$node" =~ dali ]]; then
      diffReport=diffreport.x
   else
      diffReportinPath=`command -v diffreport`
      if [ "$diffReportinPath" == "" ]; then
	 echo " ------ diffreport executable not found..."
	 echo " ------ Will use _diff_ but final report may show a FAILURE."
	 diffReport="diff -q"
      else
         diffReport="diffreport"
      fi
   fi
}

# -------------------------------------------------------------------
initEnvironment()
# -------------------------------------------------------------------
{
   begTime="$(date +%s)"
   
   diagMessage "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_"
   diagMessage "Date: $(date)"
   diagMessage "Host: $(hostname)"
   diagMessage "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_"
   diagMessage " "
   diagMessage " *** Setup $rundeck testing environment ***"

# Checkpoint initial path
   startDirectory=`pwd`

# set some directory paths
   baselineAndDataPaths

#if on DISCOVER, update git repository
   cd $workPath
   if [[ "$node" =~ borg ]]; then
      git pull
   fi

# set email address and log file names
   setLogNames

# check/set module environment
   checkSetModEnv "$moduleSet"

# set NPES based on rundeck
   setNPES

# Set EXTRA_FFLAGS 
   if [ "$fflags" == "NO" ]; then
      EXTRA_FFLAGS=""
   fi

# If we are doing baseline run then make sure of these defaults
   if [ "$baseline" == "YES" ]; then
      check="NO"
      clean="YES"
   fi

# exclude PFUNIT
   unset PFUNIT

   printInfo
}

# -------------------------------------------------------------------
setLogNames()
# -------------------------------------------------------------------
{
# log file
   makeLog=${workPath}/"make-modele-log"
   rm -f $makeLog
}

# -------------------------------------------------------------------
checkSetModEnv()
# -------------------------------------------------------------------
{
   local moduleSet=$1
   local modList

# On DISCOVER compute nodes
   if [[ "$node" =~ borg ]]; then
      diagMessage " --- Running in $pbsType mode"
      diagMessage " ------ setModules $moduleSet"
      setModules "$moduleSet"
   else
      diagMessage " ------ Not on nccs.nasa.gov systems."
      if [ "$rcFile" == "" ]; then
         diagMessage " ------ rcFile not specified. Will use $HOME/.modelErc"
         rcFile=$HOME/.modelErc
      fi
      modList=`module list > /dev/null 2>&1`
      ifcCheck=`which ifort`
      gccCheck=`which gfortran`
      mpiCheck=`which mpirun`
      if [ "$modList" == "No Modulefiles Currently Loaded." ]; then
         diagMessage " ------ No Modulefiles Currently Loaded."
      fi  
      if [ "$ifcCheck" == "" ]; then
         diagMessage " ------ Intel compiler not found."
      else
         diagMessage " ------ Intel compiler : YES"
      fi  
      if [ "$gccCheck" == "" ]; then
         diagMessage " ------ gfortran compiler not found."
      else
         diagMessage " ------ gfortran compiler : YES"
      fi  
      if [ "$mpiCheck" == "" ]; then
         diagMessage " ------ No MPI installation found."
         diagMessage " ------ WARNING: MPI compilation and/or model execution will FAIL"
      else
         diagMessage " ------ MPIRUN : "$mpiCheck
         mpiCheck=`mpif90 --version | head -1`
         diagMessage " ------ MPIF90 will use the following compiler : "$mpiCheck
      fi  

# Since we could have two compilers installed, set EXTRA_FFLAGS based on MPIcheck
      if [[ "$mpiCheck" =~ GNU ]]; then
         EXTRA_FFLAGS="-O0 -g -fbacktrace"
      else # intel
         EXTRA_FFLAGS="-O0 -g -traceback"
      fi  
      gitCheck=`which git`
      if [ "$gitCheck" == "" ]; then
         finalize 1 " --- No GIT installation found."
      fi
   fi
}

# -------------------------------------------------------------------
baselineAndDataPaths()
# -------------------------------------------------------------------
{
# May need to adjust gitBaseRepository and/or gcmSearchPath

   if [[ "$check" == "YES"  || "$baseline" == "YES" ]]; then
      if [ -z "$MODELEBASE" ]; then
         if [[ "$node" =~ borg ]]; then
            gitBaseRepository=$NOBACKUP/devel/modelE.clones/base
         else
            gitBaseRepository=$HOME/models/devel/modelE.clones/base
         fi
         diagMessage " ------ MODELEBASE was not set. Will use $gitBaseRepository"
      else
         gitBaseRepository=$MODELEBASE
      fi
   fi
   if [[ "$node" =~ borg ]]; then
      gcmSearchPath=/discover/nobackup/projects/giss/prod_input_files
   else
      gcmSearchPath=$HOME/data/modelE
   fi
   if [ "$baseline" == "YES" ]; then
      workPath=$gitBaseRepository/decks
   else
      workPath=$gitRepository/decks
   fi
}

# -------------------------------------------------------------------
setModules()
# -------------------------------------------------------------------
{
# On DISCOVER, set module environment based  on input argument. Here we also set
# the library paths used in .modelErc

    local moduleSet=$1

    source /usr/share/modules/init/bash
    module purge
    # NOTE: DEFAULT is intel11
    if [ "$moduleSet" == "intel13" ]; then
       module load comp/intel-13.0.0.079 mpi/impi-3.2.2.006
       PNETCDFHOME=/usr/local/other/pnetcdf/intel12.0.1.107_impi3.2.2.006
       NETCDFHOME=/usr/local/other/netcdf/3.6.2_intel-12.0.1.107
       ESMFHOME=/usr/local/other/esmf400rp1/intel12_impi32
       EXTRA_FFLAGS="-O0 -g -traceback"
    elif [ "$moduleSet" == "gcc47" ]; then
       module load other/comp/gcc-4.7-20120331 other/mpi/mvapich2-1.8a2/gcc-4.7-20120331
       PNETCDFHOME=/usr/local/other/esmf400rp1/gcc4.7_mvapich2-1.8
       NETCDFHOME=/usr/local/other/pnetcdf/gcc4.6_mvapich2-1.6
       ESMFHOME=/usr/local/other/esmf400rp1/gcc4.7_mvapich2-1.8
       COMPILER=gfortran
       MPIDISTR=mvapich2
       EXTRA_FFLAGS="-O0 -g -fbacktrace"
    else
      finalize 1 " --- Module set $moduleSet is not available."
    fi
    # finally load GIT module
    module load other/git-1.7.3.4
}

# -------------------------------------------------------------------
setNPES()
# -------------------------------------------------------------------
{
# In this function we determine, based on the rundeck, how many
# cpus to use. 
   if [[ "$rundeck" =~ EM20  || "$rundeck" =~ E_AR5_C12 ]]; then
      npes=( 4 )
   elif [[ "$rundeck" =~ E4Tcad || "$rundeck" =~ arobio  ]]; then
      npes=( 8 )
   elif [[ "$rundeck" =~ E4C90 ]]; then
      npes=( 6 )
   else 
      npes=( 8 )
   fi
}

# -------------------------------------------------------------------
buildAndRun() 
# -------------------------------------------------------------------
{
   local codebase
   echo " *** Processing $rundeck ... please wait ***"
   diagMessage " ------------------------------------------"
   if [ "$baseline" == "YES" ]; then
      codebase="base"
   else
      codebase=exp
      diagMessage " --- EXPERIMENTAL ---"
   fi
   buildAndRunCodeBase "$codebase"
   if [[ "$baseline" == "YES" ]]; then finalize 0; fi
}

# -------------------------------------------------------------------
buildAndRunCodeBase() 
# -------------------------------------------------------------------
{
   local codebase=$1
   local val
   local cnt
   local i
   if [ "$baseline" == "YES" ]; then

      val=${npes[0]}
      make vclean 1>> $makeLog 2>&1 

      rmFile "${rundeck}.${codebase}.exe"
      wait
      buildMpi "${rundeck}.${codebase}" 
      runModel "${rundeck}.${codebase}" "$val"
      restartRegression "${rundeck}.${codebase}" "$val"

   else

# Remove exe file (just in case) 
      rmFile "${rundeck}.${codebase}.exe"
 
      if [[ "$codebase" != "base" ]]; then
         cnt=0
         buildOption=YES
         for val in "${npes[@]}"; do
   	    if [ $cnt -gt 0 ]; then
               buildOption="NO"
     	    fi 
	    buildMpi "${rundeck}.${codebase}" "$buildOption"
	    let cnt=$cnt+1
	    runModel "${rundeck}.${codebase}" "$val"
#         saveState "${rundeck}.${codebase}" "$val"
# Run restart regression (compares 25 hr restart vs 24+1hr restart)
            if [ "$regress" == "YES" ]; then
               restartRegression "${rundeck}.${codebase}" "$val"
            fi
         done
      fi
   fi
}

# -------------------------------------------------------------------
buildSerial() 
# -------------------------------------------------------------------
{
   local deck=$1

   diagMessage " --- Building modelE SERIAL version"

# make rundeck
   makeRundeck "$deck"
   exitIfNofile "$deck.R"

   if [ "$clean" == "YES" ]; then   
      make --quiet vclean 1>> $makeLog 2>&1
   fi

# make gcm
   echo " ====== START SERIAL BUILD ====== " >> $makeLog
   diagMessage " ------ make -j gcm RUN=$deck"
   if [ "$EXTRA_FFLAGS" == "" ]; then   
      make -j gcm RUN=$deck MODELERC=$rcFile 1>> $makeLog 2>&1
   else
      make -j gcm RUN=$deck MODELERC=$rcFile EXTRA_FFLAGS="$EXTRA_FFLAGS" 1>> $makeLog 2>&1
   fi
   wait
   exitIfNofile "$deck.exe" "COMPILATION failed" 
   diagMessage " --- SERIAL compilation successful"
}

# -------------------------------------------------------------------
buildMpi() 
# -------------------------------------------------------------------
{
   local deck=$1
   local buildOpt=$2

   if [ "$baseline" == "YES" ]; then

      diagMessage " --- Building modelE baseline MPI version"

# make rundeck
      makeRundeck "$deck"
      exitIfNofile "$deck.R"

      make --quiet vclean 1>> $makeLog 2>&1
# make gcm
      echo " ====== START MPI BUILD ====== " >> $makeLog
      diagMessage " ------ make -j grucm RUN=$deck MPI=YES" 
      if [ "$EXTRA_FFLAGS" == "" ]; then   
         make -j gcm RUN=$deck MODELERC=$rcFile MPI=YES 1>> $makeLog 2>&1
      else
         make -j gcm RUN=$deck MODELERC=$rcFile MPI=YES EXTRA_FFLAGS="$EXTRA_FFLAGS" 1>> $makeLog 2>&1
      fi

   else

      if [ "$clean" == "YES" ]; then   
         make --quiet vclean 1>> $makeLog 2>&1
      fi

# Set buildOpt to avaoid recompiling
      if [ "$buildOpt" == "YES" ]; then

         diagMessage " --- Building modelE MPI version"

# make rundeck
         makeRundeck "$deck"
         exitIfNofile "$deck.R"

# make gcm
         echo " ====== START MPI BUILD ====== " >> $makeLog
         diagMessage " ------ make -j gcm RUN=$deck MPI=YES" 
         if [ "$EXTRA_FFLAGS" == "" ]; then   
            make -j gcm RUN=$deck MODELERC=$rcFile MPI=YES 1>> $makeLog 2>&1
         else
            make -j gcm RUN=$deck MODELERC=$rcFile MPI=YES EXTRA_FFLAGS="$EXTRA_FFLAGS" 1>> $makeLog 2>&1
         fi

      fi # end if buildOpt
   fi
   wait
   exitIfNofile "$deck.exe" "COMPILATION failed" 
   diagMessage " --- MPI compilation successful"

}

# -------------------------------------------------------------------
makeRundeck() 
# -------------------------------------------------------------------
{ 
   local deck=$1
   #if [ ! -e ../templates/$rundeck.R ]; then
   #   finalize 1 " --- RUNSRC=$rundeck does not exist."
   #fi
   diagMessage " ------ make rundeck RUN=$deck RUNSRC=$rundeck"
   make rundeck RUN=$deck RUNSRC=$rundeck OVERWRITE=YES MODELERC=$rcFile 1>> $makeLog 2>&1
}

# -------------------------------------------------------------------
runModel() 
# -------------------------------------------------------------------
{   
   local deck=$1
   local ncpu=$2
   if [ $ncpu -eq 0 ]; then
      diagMessage " --- Run $deck (serial mode)"
      make -j setup RUN=$deck MODELERC=$rcFile 1>> $makeLog 2>&1
      MODELERC=$rcFile ../exec/runE $deck -np 1 -cold-restart 1>> $makeLog 2>&1
   else
       diagMessage " --- Run $deck (MPI mode,  NPES=$ncpu)"
       make -j setup RUN=$deck MODELERC=$rcFile MPI=YES 1>> $makeLog 2>&1
       MODELERC=$rcFile ../exec/runE $deck -np $ncpu -cold-restart 1>> $makeLog 2>&1
   fi
   checkRun "$deck"
}

# -------------------------------------------------------------------
checkRun() 
# -------------------------------------------------------------------
{   
   local deck=$1
   if [ -e $deck/run_status ]; then
      rcCode=`head -1 $deck/run_status`
      rcMess=`tail -1 $deck/run_status`
      diagMessage " ------ run_status : rc=$rcCode; $rcMess"
      if [ $rcCode -ne 13 ]; then
         finalize 1 " ------ Problem encountered while running $deck";
      fi
   else
      finalize 1 " ------ $deck/run_status does not exist";
   fi 
}

# -------------------------------------------------------------------
restartRegression() 
# -------------------------------------------------------------------
{ 
   local deck=$1
   local ncpu=$2
   local regDeck=${deck}.reg
   diagMessage " --- Regression (1day) runs"

# Make a new rundeck - it will have the .reg suffix
   makeRundeck "$regDeck"

#   if [ "$baseline" == "YES" ]; then

# Edit rundeck to run 1 day
# Note that 48 = # time steps in one day
      editRundeck "$regDeck.R" 48 2 0
      diagMessage " --- Run model one day"
      runModel "$regDeck" "$ncpu"   
  
#   else

#      if [ "$regress" == "YES" ]; then
# Edit rundeck to run 1 hr past day 2 (any assumptions here?) and
# create a checkpoint (fort.1) at day 2.
#         editRundeck "$regDeck.R" 48 2 1
# Run model from t0 -> t2 (running one day)
#         runModel "$regDeck" "$ncpu"
# Remove run_status - this is only a problem on discover (why?)
#         rm -f $regDeck/run_status
# Save t2 (fort.2) into t1 (fort.1) , and keep a copy of t2 for
# comparison with t2'
#         saveState "$regDeck" "$ncpu" "YES"
#         cd $regDeck
# Run model from from t1 -> t2'
#         diagMessage " --------- ./$regDeck -np $ncpu"
#         ./$regDeck -np "$ncpu" 1>> $makeLog 2>&1
#      fi

#   fi
   cd $workPath

}

# -------------------------------------------------------------------
editRundeck()
# -------------------------------------------------------------------
{
# Adapted from R. Ruedy's script add_regr_test_params
# This function enables regression testing of restarts. It modifies
# the RUNDECK by adding an NDISK line right before the &&END_PARAMETERS 
# and an end time near the end of the RUNDECK.

   local deck=$1
   local ndisk=$2
   local datee=$3
   local houre=$4

   ndisk_line="ndisk=$ndisk"
   end_hour_line=" DATEE=$datee, HOURE=$houre,"
   eof1='&&END_PARAMETERS'
   eof2='/'

   a=$( grep -n ${eof1}  $deck | head -1 ) ; n1=${a%%:*}
   a=$( grep -n ^.${eof2} $deck | head -1 ) ; n2=${a%%:*}

   cp ${deck} templ
   head -$(( n1-1 )) templ                   > ${deck}
   echo "${ndisk_line}"                      >> ${deck}
   tail +${n1} templ | head -$(( n2 - n1 ))  >> ${deck}
   echo "${end_hour_line}"                   >> ${deck}
   echo " ISTART=2, ${end_hour_line}"         >> ${deck}
   tail +${n2} templ                         >> ${deck}
}

# -------------------------------------------------------------------
compareRuns() 
# -------------------------------------------------------------------
{
   local n
   local cnt
   local rc
   local base=$gitBaseRepository/decks/$rundeck.base
   local exp=$workPath/$rundeck.exp
   local baseDay=$gitBaseRepository/decks/$rundeck.base.reg
   local exp=$workPath/$rundeck.exp
   local expDay=$workPath/$rundeck.exp.reg
   local expMpi=$workPath/$rundeck.exp
   local expMpiDay=$workPath/$rundeck.exp.reg
   declare -a id
   id=( )

   if [ ! -L $base ]; then
      finalize 1 " ------ $base does not link to an output directory";
   fi

   if [ "$compile" == "YES" ]; then
      return
   fi

   echo " *** Results ***"
   echo " ------------------------------------------"
   cnt=0
   for n in "${npes[@]}"; do
      id[${#id[*]}]="NPES${n}"
      let cnt=$cnt+1
   done

#   echo " --- DIFF results: Strong reproducibility"   
#   if [ "$runSerial" == "YES" ]; then
#      for (( n=0 ;  n<${#npes[*]}; n++ )); do    
#         rc=`$diffReport $exp/fort.2.nc $expMpi/fort.2.nc.${id[$n]} > .sz`
#         sz=`ls -la .sz | awk '{print $5}'`
#         checkStatus $sz " ------ SERIAL vs ${id[$n]}: "
#         rm -f .sz
#      done
#   fi

#   let cnt=${#npes[*]}-1
#   for (( n=0 ;  n<$cnt; n++ )); do
#      let n1=$n+1
#      rc=`$diffReport $expMpi/fort.2.nc.${id[0]} $expMpi/fort.2.nc.${id[$n1]} > .sz`
#      sz=`ls -la .sz | awk '{print $5}'`
#      checkStatus $sz " ------ ${id[0]} vs ${id[$n1]}: "
#      rm -f .sz
#   done

   if [ "$check" == "YES" ]; then
      echo " --- DIFF results: Has model changed?"   
         exitIfNofile "$base/fort.2.nc"
         exitIfNofile "$exp/fort.2.nc"
         rc=`$diffReport $base/fort.2.nc $exp/fort.2.nc > .sz`
         sz=`ls -la .sz | awk '{print $5}'`
         checkStatus $sz " ------ BASELINE vs EXP 1HR: "
         rm -f .sz

         exitIfNofile "$baseDay/fort.2.nc"
         exitIfNofile "$expDay/fort.2.nc"
         rc=`$diffReport $baseDay/fort.2.nc $expDay/fort.2.nc > .sz`
         sz=`ls -la .sz | awk '{print $5}'`
         checkStatus $sz " ------ BASELINE vs EXP 1DY: "
         rm -f .sz

#      let cnt=${#npes[*]}
#      for (( n=0 ;  n<$cnt; n++ )); do
#         rc=`$diffReport $base/fort.2.nc $expMpi/fort.2.nc.${id[$n]} > .sz`
#         sz=`ls -la .sz | awk '{print $5}'`
#         checkStatus $sz " ------ BASELINE vs EXP (${id[$n]}): "
#         rm -f .sz
#      done
   fi
   unset id
   echo " ------------------------------------------"
}

# -------------------------------------------------------------------
compareRegressionRuns() 
# -------------------------------------------------------------------
{
   local n
   local cnt
   local rc
   local out2=$rundeck.exp.reg
   declare -a id
   id=( )

   cnt=0
   for n in "${npes[@]}"; do
      id[${#id[*]}]="NPES${n}"
      let cnt=$cnt+1
   done

   echo " --- DIFF results - Restart regression"   

   let cnt=${#npes[*]}
   for (( n=0 ;  n<$cnt; n++ )); do
      echo $diffReport $out2/restart.1dy.${id[$n]} $out2/fort.2.nc
      rc=`$diffReport $out2/restart.1dy.${id[$n]} $out2/fort.2.nc > .sz`
      sz=`ls -la .sz | awk '{print $5}'`
      checkStatus $sz " ------ restart regression for ${id[$n]} : "
      rm -f .sz
   done
   unset id
   echo " ------------------------------------------"
}

# -------------------------------------------------------------------
saveState() 
# -------------------------------------------------------------------
{
   local deck=$1
   local n=$2
   local saveReg=$3
   local id="NPES${n}"
   diagMessage " ------ saveState $deck $id"
# Save checkpoints
   if [ "$saveReg" == "YES" ]; then
      cp $deck/fort.1.nc $deck/restart.1dy.$id
      cp $deck/fort.2.nc $deck/fort.1.nc
   else
      cp $deck/fort.2.nc $deck/fort.2.nc.$id
   fi
}

# -------------------------------------------------------------------
finalize()
# -------------------------------------------------------------------
# Mail report and clean up
{
   local exitCode=$1
   local exitMessage=$2
   if [ $exitCode -eq 1 ]; then
      diagMessage " @@@ $rundeck run: FAILURE @@@"   
      diagMessage "$exitMessage"
      exit 1
   fi
   diagMessage " *** $rundeck run is complete. ***"

   endTime="$(date +%s)"
   elapsed="$(expr $endTime - $begTime)"
   remainder="$(expr $elapsed % 3600)"
   hours="$(expr $(expr $elapsed - $remainder) / 3600)"
   seconds="$(expr $remainder % 60)"
   minutes="$(expr $(expr $remainder - $seconds) / 60)"

   diagMessage " "
   diagMessage "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_"
   echo "Total time: $hours:$minutes:$seconds"
   diagMessage "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_"
   cd $startDirectory
   exit 0
}

# -------------------------------------------------------------------
checkStatus() 
# -------------------------------------------------------------------
{
   local sz=$1
   local msg=$2

   if [[ $sz -eq 0 ]]; then
      echo "$msg""0-diff"
      return
   fi
   echo "$msg"" NON-ZERO diff <------ FAILURE"
   finalize 1 "checkStatus gives non-zero diffs"
}

# -------------------------------------------------------------------
exitIfNofile() 
# -------------------------------------------------------------------
{ 
   local fullpath=$1
   local msg=$2
   filename="${fullpath##*/}" 
   dir="${fullpath:0:${#fullpath} - ${#filename}}" 
   if [ "$dir" == "" ]; then dir="."; fi
   result=`find $dir -name $filename`
   if [ "$result" == "" ]; then
      if [ "$msg" != "" ]; then
         diagMessage " ------ $msg"
         finalize 1 " ------ Please see $makeLog"
      else
         finalize 1 " ------ $fullpath does not exist."
      fi
   fi
}

# -------------------------------------------------------------------
rmFile() 
# -------------------------------------------------------------------
{ 
   local file=$1
   find . -name $file -exec rm -f {} \; > /dev/null 2>&1
}

# -------------------------------------------------------------------
printInfo()
# -------------------------------------------------------------------
{   
   diagMessage " *** Options ***"
   diagMessage "    EXPERIMENTAL GIT repository: $gitRepository"
   if [[ "$check" == "YES"  || "$baseline" == "YES" ]]; then
      diagMessage "    create baseline               = $baseline"
      diagMessage "    BASELINE GIT repository: $gitBaseRepository"
   fi
   if [[ "$pbsType" =~ BATCH ]]; then
      diagMessage "    moduleEnv                     = $moduleSet"
   fi
   diagMessage "    modelErc file                 = $rcFile"
   diagMessage "    runDeck                       = $rundeck"
   diagMessage "    make clean                    = $clean"
   diagMessage "    restart regression            = $regress"
   diagMessage "    check against base repository = $check"
   if [ "$fflags" != "" ]; then
      diagMessage "    fortran flags                 = $EXTRA_FFLAGS"
   fi
}

# -------------------------------------------------------------------
diagMessage()
# -------------------------------------------------------------------
{
   local text=$1
   if [ "$verbose" == "YES" ]; then
      echo "$text"
   fi
}

# -------------------------------------------------------------------
parseOptions()
# -------------------------------------------------------------------
{
# if these are defined thru ENV then set them; but command line will override them

   if [ ! -z "$MODELERC" ]; then rcFile=$MODELERC; fi
   if [ ! -z "$MODELEDECK" ]; then rundeck=$MODELEDECK; fi

# Loop through argument options

   for (( i = 0 ; i < ${#names[@]} ; i++ ))
      do

       if [[ "${names[$i]}" =~ --modEnv ]]; then 
           moduleSet=`echo ${names[$i]} | sed 's/[-a-zA-Z0-9]*=//'`
       elif [[ "${names[$i]}" =~ -e ]]; then 
           let i=i+1
	   moduleSet=${names[$i]}
	   
       elif [[ "${names[$i]}" =~ --deck ]]; then 
	   rundeck=`echo ${names[$i]} | sed 's/[-a-zA-Z0-9]*=//'`
       elif [[ "${names[$i]}" =~ -d ]]; then 
           let i=i+1
	   rundeck=${names[$i]}
	   
       elif [[ "${names[$i]}" =~ --modelErc ]]; then 
	   rcFile=`echo ${names[$i]} | sed 's/[-a-zA-Z0-9]*=//'`
       elif [[ "${names[$i]}" =~ -m ]]; then 
           let i=i+1
	   rcFile=${names[$i]}
           
       elif [[ "${names[$i]}" =~ --fflags ]]; then 
	   fflags=`echo ${names[$i]} | sed 's/[-a-zA-Z0-9]*=//'`
       elif [[ "${names[$i]}" =~ -f ]]; then 
           let i=i+1
	   fflags=${names[$i]}
	   
       elif [[ "${names[$i]}" =~ --verbose ]]; then 
	   verbose=`echo ${names[$i]} | sed 's/[-a-zA-Z0-9]*=//'`
       elif [[ "${names[$i]}" =~ -v ]]; then 
           let i=i+1
	   verbose=${names[$i]}
	   
       elif [[ "${names[$i]}" =~ --clean ]]; then 
	   clean=`echo ${names[$i]} | sed 's/[-a-zA-Z0-9]*=//'`
       elif [[ "${names[$i]}" =~ -c ]]; then 
           let i=i+1
	   clean=${names[$i]}
	   
       elif [[ "${names[$i]}" =~ --check ]]; then 
	   check=`echo ${names[$i]} | sed 's/[-a-zA-Z0-9]*=//'`
       elif [[ "${names[$i]}" =~ -x ]]; then 
           let i=i+1
	   check=${names[$i]}
	   
       elif [[ "${names[$i]}" =~ --regress ]]; then 
   	   regress=`echo ${names[$i]} | sed 's/[-a-zA-Z0-9]*=//'`
       elif [[ "${names[$i]}" =~ -r ]]; then 
           let i=i+1
	   regress=${names[$i]}

       elif [[ "${names[$i]}" =~ --baseline ]]; then 
   	   baseline=`echo ${names[$i]} | sed 's/[-a-zA-Z0-9]*=//'`
       elif [[ "${names[$i]}" =~ -b ]]; then 
           let i=i+1
	   baseline=${names[$i]}

       else
           echo " ### Unknown argument ${names[$i]}"; exit 1
       fi
   done # end do loop
}

# Command line parser:
# -------------------------------------------------------------------
getopt()
# -------------------------------------------------------------------
{
  var=""
  wantarg=0
  for (( i=1; i<=$#; i+=1 )); do
    lastvar=$var
    var=${!i}
    if [ "$var" = "" ]; then 
        continue 
    fi
    echo \ $var | grep -q -- '='
    if [ $? -eq 0 ]; then
      ## -*param=value
      var=$(echo \ $var | sed -r s/'^[ ]*-*'/''/)
      myvar=${var%=*}
      myval=${var#*=}
      eval "${myvar}"="'$myval'"
    else
      echo \ $var | grep -E -q -- '^[ ]*-'
      if [ $? -eq 0 ]; then
        ## -*param$
        var=$(echo \ $var | sed -r s/'^[ ]*-*'/''/)
        eval "${var}"=1
        wantarg=1
      else
        echo \ $var | grep -E -- '^[ ]*-'
        if [ $? -eq 0 ]; then
          # the current one has a dash, so cannot be
          # the argument to the last parameter
          wantarg=0
        fi
        if [ $wantarg -eq 1 ]; then
          # parameter argument
          val=$var
          var=$lastvar
          eval "${var}"="'${val}'"
          wantarg=0
        else
          # parameter
          if [ "${!var}" = "" ]; then
            eval "${var}"=1
          fi
          wantarg=0
        fi
      fi
    fi
  done
}

# Trap not-normal exit signals: 1/HUP, 2/INT, 3/QUIT, 15/TERM
# @see catch_sig()
trap catch_sig 1 2 3 15

# -------------------------------------------------------------------
function cleanexit() 
#  Wrapper around 'exit' to cleanup on exit.
#  @param $1 integer  Exit status.  If $1 not defined, exit status of global
#+                    variable 'EXIT_STATUS' is used.  If neither $1 or
#+                    'EXIT_STATUS' defined, exit with status 0 (success).
# -------------------------------------------------------------------
{
    echo "Exiting with ${1:-${EXIT_STATUS:-0}}"
    exit ${1:-${EXIT_STATUS:-0}}
} 

# -------------------------------------------------------------------
function catch_sig() 
# Catch signal trap.
# Trap not-normal exit signals: 1/HUP, 2/INT, 3/QUIT, 15/TERM
# @NOTE1: Non-trapped signals are 0/EXIT, 9/KILL.
# -------------------------------------------------------------------
{
    local exit_status=$?
    echo "$0 terminated abnormally..."
    cleanexit $exit_status
} 
