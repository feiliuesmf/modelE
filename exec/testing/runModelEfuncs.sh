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
     environment variable. If specified, the results will be compared against a set of
     baseline results. (Baseline results must be precomputed.)

     To create a modelE baseline just set the MODELEBASE environment variable and the
     corresponding command line flag.
    
     On DISCOVER it is assumed that the script is executed within a PBS interactive
     session. On non-DISCOVER systems a fortran compiler (gfortran or Intel) and
     MPI must be available.    
 

SYNOPSIS

     $0 <...arguments...>  [..options...]

REQUIREMENTS

     MODELEGIT      : environment variable defining git local repository path
                      e.g. export MODELEGIT=<path to git repository>

ARGUMENTS

     -h | --help     : this help screen

OPTIONS

     MODELEBASE      : environment variable defining a baseline repository path
                      e.g. export MODELEBASE=<path to git repository>

     -b | --baseline : create modelE baseline (default: NO) 
     -c | --clean    : run "make clean" (default: YES) 
     -d | --deckname : use a different name for rundeck (default: same as template name)
     -f | --fflags   : use debug FFLAGS to compile (default: YES)
     -m | --modelErc : specified location of .modelErc (default: $HOME/.modelErc)
                       OR MODELERC environment variable
     -o | --compile  : compile only (default: NO)
     -r | --runtype  : runtype is 1HR or 1DY or ALL (default: ALL)
     -l | --log      : view build/model output (default: NO)
     -t | --template : rundeck template name (default: nonProduction_E_AR5_C12)
     -x | --check    : check run results against baseline repository (default: NO)
                       If yes, MODELEBASE environement variable must be set
     -v | --verbose  : verbose (default: YES)

     Additionally, on DISCOVER:

     -e | --modEnv   : module environment (default: intel)
                             intel : comp/intel-13.0.0.079 
                                     mpi/impi-3.2.2.006
                             gcc   : other/comp/gcc-4.7-20120331
				     other/mpi/mvapich2-1.8a2/gcc-4.7-20120331
EXAMPLE

     This screen:
        runModelE.sh -h 
     Create baseline (with MODELEBASE set):
        runModelE.sh
     Do not check against baseline (just runs within MODELEGIT)
        runModelE.sh -x NO
     Use gcc modules, E4TcadF40 rundeck (with shortname=cad) and the specified modelErc file
        runModelE.sh --modEnv=gcc -t E4TcadF40 -d cad --modelErc=$HOME/.modelErc.G47

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

# .modelErc file path - the default file is automatically created here
   rcFile=$HOME/.modelErc

# Rundeck name
   template="nonProduction_E_AR5_C12"

# Default rundeck name is the same as the template name
   deckname="nonProduction_E_AR5_C12"

# module environment in .modelErc - default is to use Intel v13.x
   moduleSet="intel"
   # used in .modelErc:
   COMPILER=intel
   MPIDISTR=intel
   MPIDIR=

# Debug compilation flags
   fflags="YES"

# Create a set of modelEbaseline restarts to be compared against
   baseline="NO"

# Run make clean
   clean="YES"

# Compile only
   compile="NO"

# Run type ALL=1hr+1dy+regression
   runtype="ALL"

# Print out run information
   verbose="YES"

# Print out run information
   viewlog="NO"

# Check rundeck changes against those in MODELEBASE repository
   check="NO"

# "diff report" for nc files based on GISS utility. On DISCOVER it is 
# installed and ready to go but user must build it on other systems (and
# make sure it is in the PATH).
   if [[ "$node" =~ borg || "$node" =~ discover || "$node" =~ dali ]]; then
      diffReport=diffreport.x
   else
      diffReportinPath=`command -v diffreport.x`
      if [ "$diffReportinPath" == "" ]; then
	 echo " ------ diffreport executable not found..."
	 echo " ------ Will use _diff_ but final report may show a FAILURE."
	 diffReport="diff -q"
      else
         diffReport="diffreport.x"
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
   diagMessage " *** Setup $template testing environment ***"

# Checkpoint initial path
   startDirectory=`pwd`

# set some directory paths
   baselineAndDataPaths

#if on DISCOVER, update git repository (disabling for now)
   cd $workPath
   #if [[ "$node" =~ borg && "$baseline" == "NO" ]]; then
   #   git pull > /dev/null 2>&1
   #fi

# set email address and log file names
   setLogNames

# check/set module environment
   checkSetModEnv "$moduleSet"

# Set EXTRA_FFLAGS 
   if [ "$fflags" == "NO" ]; then
      EXTRA_FFLAGS=""
   fi

# For convenience define a glocal variable named rundeck
   rundeck=$template
   if [ "$deckname" != "$template" ]; then
      rundeck=$deckname
   fi

# set NPES based on rundeck
   setNPES

# If we are doing baseline run then ...
   if [ "$baseline" == "YES" ]; then
      check="NO"
#      clean="YES"
   fi

# exclude PFUNIT
   unset PFUNIT

   export MODELERC=$rcFile

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
# if we are checking against baseline or doing a baseline run then we need
# a baseline repository...first check if it's valid, else use some defaults
   if [[ "$check" == "YES" || "$baseline" == "YES" ]]; then
      if [ -z "$MODELEBASE" ]; then
         if [[ "$node" =~ borg ]]; then
           # WARNING: This setting is not meant to work for everyone
           # Use MODELEBASE instead!
            gitBaseRepository=$NOBACKUP/devel/modelE.clones/base
         else
           # WARNING: If working on your own system then change this or
           # use MODELEBASE!
            gitBaseRepository=$HOME/models/devel/modelE.clones/base
         fi
         diagMessage " ------ MODELEBASE was not set. Will use $gitBaseRepository"
      else
         gitBaseRepository=$MODELEBASE
         if [ ! -d "$gitBaseRepository" ]; then
            finalize 1 " --- $gitBaseRepository does not exist"
         fi
      fi
   fi
   if [[ "$node" =~ borg ]]; then
      gcmSearchPath=/discover/nobackup/projects/giss/prod_input_files
   else
      # WARNING: If working on your own system then change this:
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
    # NOTE: DEFAULT is intel
    if [ "$moduleSet" == "intel" ]; then
       module load comp/intel-13.0.0.079 mpi/impi-3.2.2.006
       PNETCDFHOME=/usr/local/other/pnetcdf/intel12.0.1.107_impi3.2.2.006
       NETCDFHOME=/usr/local/other/netcdf/3.6.2_intel-12.0.1.107
       ESMFHOME=/usr/local/other/esmf400rp1/intel12_impi32
       EXTRA_FFLAGS="-O0 -g -traceback"
    elif [ "$moduleSet" == "gcc" ]; then
       module load other/comp/gcc-4.7.1 other/mpi/mvapich2-1.9a2/gcc-4.7.1
       PNETCDFHOME=/usr/local/other/pnetcdf/gcc-4.7.1_mvapich2-1.9a2
       NETCDFHOME=/usr/local/other/netcdf/3.6.2_gcc4.6
       ESMFHOME=/usr/local/other/esmf400rp1/gcc4.7_mvapich2-1.8
       COMPILER=gfortran
       MPIDISTR=mvapich2
       MPIDIR=/usr/local/other/mvapich2/1.9a2/gcc-4.7.1
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
   if [ "$npes" == "" ]; then
     if [[ "$template" =~ cadC12  || "$template" =~ AR5_C12 ]]; then
       useNPES=( 2 )
     elif [[ "$template" =~ EM20 || "$template" =~ E1oM20 ]]; then
       useNPES=( 4 )
     elif [[ "$template" =~ cadiF40 || "$template" =~ arobio  ]]; then
       useNPES=( 8 )
     elif [[ "$template" =~ E4C90 ]]; then
       useNPES=( 24 )
     elif [[ "$template" =~ tomas || "$template" =~ campiF40 || "$template" =~ ampF40 ]]; then
       useNPES=( 44 )
     else 
       useNPES=( 12 )
     fi
   else
      useNPES=$npes
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
      codebase=dev
      diagMessage " --- DEVELOPMENT ---"
   fi
   buildAndRunCodeBase "$codebase"
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

      val=${useNPES[0]}
#      if [ "$clean" == "YES" ]; then   
#         make vclean 1>> $makeLog 2>&1 
#      fi

      rmFile "${rundeck}.${codebase}.exe"
      wait
      buildMpi "${rundeck}.${codebase}" 
      if [ "$compile" == "YES" ]; then
        return 
      fi
      runModel "${rundeck}.${codebase}" "$val"
      saveState "${rundeck}.${codebase}" "$val" "HRRUN"

      if [ "$runtype" == "ALL" ]; then
        restartRegression "${rundeck}.${codebase}" "$val"
      elif [ "$runtype" == "1DY" ]; then
        run1Day "${rundeck}.${codebase}" "$val"
      elif [ "$runtype" == "1HR" ]; then
        return 
      fi
   else

# Remove exe file (just in case) 
      rmFile "${rundeck}.${codebase}.exe"
 
      if [[ "$codebase" != "base" ]]; then
         cnt=0
         buildOption=YES
         for val in "${useNPES[@]}"; do
   	    if [ $cnt -gt 0 ]; then
               buildOption="NO"
     	    fi 
	    buildMpi "${rundeck}.${codebase}" "$buildOption"
            if [ "$compile" == "YES" ]; then
              return 
            fi
  	    let cnt=$cnt+1

	    runModel "${rundeck}.${codebase}" "$val"
            saveState "${rundeck}.${codebase}" "$val" "HRRUN"

            if [ "$runtype" == "ALL" ]; then
              restartRegression "${rundeck}.${codebase}" "$val"
            elif [ "$runtype" == "1DY" ]; then
              run1Day "${rundeck}.${codebase}" "$val"
            elif [ "$runtype" == "1HR" ]; then
              return 
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

   if [ "$clean" == "YES" ]; then   
      make --quiet vclean 1>> $makeLog 2>&1
   fi

# make gcm
   echo " ====== START SERIAL BUILD ====== " >> $makeLog
   diagMessage " ------ make -j gcm RUN=$deck"
   if [ "$EXTRA_FFLAGS" == "" ]; then   
      make -j gcm RUN=$deck 1>> $makeLog 2>&1
   else
      make -j gcm RUN=$deck EXTRA_FFLAGS="$EXTRA_FFLAGS" 1>> $makeLog 2>&1
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

      if [ "$clean" == "YES" ]; then   
	 make --quiet vclean 1>> $makeLog 2>&1
      fi

# make gcm
      echo " ====== START MPI BUILD ====== " >> $makeLog
      diagMessage " ------ make -j gcm RUN=$deck MPI=YES" 
      if [ "$EXTRA_FFLAGS" == "" ]; then   
         if [ "$viewlog" == "YES" ]; then   
           xterm -e make -j gcm RUN=$deck MPI=YES | tee $makeLog
         else
           make -j gcm RUN=$deck MPI=YES 1>> $makeLog 2>&1
         fi
      else
         make -j gcm RUN=$deck MPI=YES EXTRA_FFLAGS="$EXTRA_FFLAGS" 1>> $makeLog 2>&1
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
            if [ "$viewlog" == "YES" ]; then   
	      xterm -e make -j gcm RUN=$deck MPI=YES | tee $makeLog
            else
              make -j gcm RUN=$deck MPI=YES 1>> $makeLog 2>&1
            fi
         else
            make -j gcm RUN=$deck MPI=YES EXTRA_FFLAGS="$EXTRA_FFLAGS" 1>> $makeLog 2>&1
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
   #if [ ! -e ../templates/$template.R ]; then
   #   finalize 1 " --- RUNSRC=$rundeck does not exist."
   #fi
   diagMessage " ------ make rundeck RUN=$deck RUNSRC=$template"
   make rundeck RUN=$deck RUNSRC=$template OVERWRITE=YES 1>> $makeLog 2>&1
   #exitIfNofile "$template.R"
}

# -------------------------------------------------------------------
runModel() 
# -------------------------------------------------------------------
{   
   local deck=$1
   local ncpu=$2

   if [ $ncpu -eq 0 ]; then
      diagMessage " --- Run $deck (serial mode)"
      make -j setup RUN=$deck 1>> $makeLog 2>&1
      ../exec/runE $deck -np 1 -cold-restart 1>> $makeLog 2>&1
   else
      diagMessage " --- Run $deck (MPI mode,  NPES=$ncpu)"
      if [ "$viewlog" == "YES" ]; then   
        xterm -e make -j setup RUN=$deck MPI=YES | tee $makeLog
        xterm -e ../exec/runE $deck -np $ncpu -cold-restart | tee $makeLog
      else
        make -j setup RUN=$deck MPI=YES 1>> $makeLog 2>&1
        ../exec/runE $deck -np $ncpu -cold-restart 1>> $makeLog 2>&1
      fi
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
run1Day() 
# -------------------------------------------------------------------
{ 
   local deck=$1
   local ncpu=$2
   diagMessage " --- 1 Day run"

# Edit rundeck to run 1 day
   editRundeck "${deck}.R" 48 2 0
# Run model (running one day)
   runModel "$deck" "$ncpu"
# Remove run_status - this is only a problem on discover (why?)
   rm -f $deck/run_status
   saveState "$deck" "$ncpu" "DAYRUN"

   cd $workPath
}

# -------------------------------------------------------------------
restartRegression() 
# -------------------------------------------------------------------
{ 
   local deck=$1
   local ncpu=$2
   diagMessage " --- Regression run"

# Edit rundeck to run 1 hr past day 2 (any assumptions here?) and
# create a checkpoint (fort.1) at day 2.
   editRundeck "${deck}.R" 48 2 1
# Run model from t0 -> t2 (running one day+1hr)
   runModel "$deck" "$ncpu"
# Remove run_status - this is only a problem on discover (why?)
   rm -f $deck/run_status
# Save t2 (fort.2) into t1 (fort.1) , and keep a copy of t2 for
# comparison with t2'
   saveState "$deck" "$ncpu" "DAYRUN"
   cd $deck
# Run model from from t1 -> t2 (running one hr continuation)
   diagMessage " --------- ./$deck -np $ncpu"
   ./$deck -np "$ncpu" 1>> $makeLog 2>&1
   saveState "$deck" "$ncpu" "RESTART"

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
   a=$( grep -n ^${eof2} $deck | head -1 ) ; n2=${a%%:*}

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
   local dev=$workPath/$rundeck.dev
   local baseDay=$gitBaseRepository/decks/$rundeck.base.reg
   local devDay=$workPath/$rundeck.dev.reg
   local devMpi=$workPath/$rundeck.dev
   local devMpiDay=$workPath/$rundeck.dev.reg
   declare -a id
   id=( )

   if [ "$check" == "YES" ]; then
      if [ ! -d $base ]; then
         finalize 1 " ------ $base does not link to an output directory";
      fi

      echo " *** Results ***"
      echo " ------------------------------------------"
      cnt=0
      for n in "${useNPES[@]}"; do
         id[${#id[*]}]="NPES${n}"
         if [[ "$runtype" == "1HR"  || "$runtype" == "ALL" ]]; then
           exitIfNofile "$base/$rundeck.base.1hr.$id"
           exitIfNofile "$dev/$rundeck.dev.1hr.$id"
         fi
         if [[ "$runtype" == "1DY"  || "$runtype" == "ALL" ]]; then
           exitIfNofile "$base/$rundeck.base.1dy.$id"
           exitIfNofile "$dev/$rundeck.dev.1dy.$id"
         fi
         if [ "$runtype" == "ALL" ]; then
           exitIfNofile "$base/$rundeck.base.restart.$id"
           exitIfNofile "$dev/$rundeck.dev.restart.$id"
         fi
         let cnt=$cnt+1
      done

      echo " --- DIFF results: Has model changed?"   
      for n in "${useNPES[@]}"; do
         id[${#id[*]}]="NPES${n}"
         if [[ "$runtype" == "1HR"  || "$runtype" == "ALL" ]]; then
           rc=`$diffReport $base/$rundeck.base.1hr.$id  $dev/$rundeck.dev.1hr.$id > .sz`
           sz=`wc -c .sz | awk '{print $1}'`
           checkStatus $sz " ------ BASELINE vs DEV 1HR ($id): "
           rm -f .sz
         fi
      
         if [[ "$runtype" == "1DY" || "$runtype" == "ALL" ]]; then
           rc=`$diffReport $base/$rundeck.base.1dy.$id  $dev/$rundeck.dev.1dy.$id > .sz`
           sz=`wc -c .sz | awk '{print $1}'`
           checkStatus $sz " ------ BASELINE vs DEV 1DY ($id): "
           rm -f .sz
         fi
         
         if [ "$runtype" == "ALL" ]; then
           rc=`$diffReport $base/$rundeck.base.restart.$id  $dev/$rundeck.dev.restart.$id > .sz`
           sz=`wc -c .sz | awk '{print $1}'`
           checkStatus $sz " ------ BASELINE vs DEV RESTART ($id): "
           rm -f .sz
         fi
      done
   fi
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

   diagMessage " ------ save state for $deck : id=$id ($saveReg)"

# Save checkpoints
   if [ "$saveReg" == "RESTART" ]; then
      # Note this is done in the $deck directory
      cp fort.2.nc $deck.restart.$id
   elif [ "$saveReg" == "DAYRUN" ]; then
      cp $deck/fort.2.nc $deck/$deck.1dy.$id
      cp $deck/fort.2.nc $deck/fort.1.nc
   else 
      cp $deck/fort.2.nc $deck/$deck.1hr.$id
   fi
}

# -------------------------------------------------------------------
finalize()
# -------------------------------------------------------------------
# Mail report and clean up
{
   local exitCode=$1
   local exitMessage=$2

   local base=$gitBaseRepository/decks/$rundeck.base
   local dev=$workPath/$rundeck.dev

   rm -f $base/fort.1.nc $base/fort.2.nc
   rm -f $dev/fort.1.nc $dev/fort.2.nc

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
   diagMessage "    DEVELOPMENT GIT repository: $gitRepository"
   if [[ "$check" == "YES"  || "$baseline" == "YES" ]]; then
      diagMessage "    BASELINE GIT repository: $gitBaseRepository"
      diagMessage "    create baseline               = $baseline"
   fi
   diagMessage "    modelErc file                 = $rcFile"
   diagMessage "    rundeck template              = $template"
   if [ "$deckname" != "TEMPLATE" ]; then
      diagMessage "    rundeck name                  = $rundeck"
   fi
   diagMessage "    make clean                    = $clean"
   diagMessage "    compile only                  = $compile"
   if [[ "$compile" == "NO" ]]; then
      diagMessage "    run type                      = $runtype"
      if [[ "$baseline" == "NO" ]]; then
         diagMessage "    check against base repository = $check"
      fi
   fi
   if [[ "$pbsType" =~ BATCH ]]; then
      diagMessage "    moduleEnv                     = $moduleSet"
   fi
   if [ "$fflags" == "YES" ]; then
      diagMessage "    fortran flags                 = $EXTRA_FFLAGS"
   fi
   diagMessage "    NPES                          = $useNPES"
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

# Loop through argument options

   for (( i = 0 ; i < ${#names[@]} ; i++ ))
      do

       if [[ "${names[$i]}" =~ --modEnv ]]; then 
           moduleSet=`echo ${names[$i]} | sed 's/[-a-zA-Z0-9]*=//'`
       elif [[ "${names[$i]}" =~ -e ]]; then 
           let i=i+1
	   moduleSet=${names[$i]}
	   
       elif [[ "${names[$i]}" =~ --template ]]; then 
	   template=`echo ${names[$i]} | sed 's/[-a-zA-Z0-9]*=//'`
       elif [[ "${names[$i]}" =~ -t ]]; then 
           let i=i+1
	   template=${names[$i]}
	   
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
	   
       elif [[ "${names[$i]}" =~ --npes ]]; then 
	   npes=`echo ${names[$i]} | sed 's/[-a-zA-Z0-9]*=//'`
       elif [[ "${names[$i]}" =~ -n ]]; then 
           let i=i+1
	   npes=${names[$i]}
	   
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
	   
       elif [[ "${names[$i]}" =~ --baseline ]]; then 
   	   baseline=`echo ${names[$i]} | sed 's/[-a-zA-Z0-9]*=//'`
       elif [[ "${names[$i]}" =~ -b ]]; then 
           let i=i+1
	   baseline=${names[$i]}

       elif [[ "${names[$i]}" =~ --compile ]]; then 
   	   compile=`echo ${names[$i]} | sed 's/[-a-zA-Z0-9]*=//'`
       elif [[ "${names[$i]}" =~ -o ]]; then 
           let i=i+1
	   compile=${names[$i]}

       elif [[ "${names[$i]}" =~ --runtype ]]; then 
   	   runtype=`echo ${names[$i]} | sed 's/[-a-zA-Z0-9]*=//'`
       elif [[ "${names[$i]}" =~ -r ]]; then 
           let i=i+1
	   runtype=${names[$i]}

       elif [[ "${names[$i]}" =~ --deckname ]]; then 
   	   deckname=`echo ${names[$i]} | sed 's/[-a-zA-Z0-9]*=//'`
       elif [[ "${names[$i]}" =~ -d ]]; then 
           let i=i+1
	   deckname=${names[$i]}

       elif [[ "${names[$i]}" =~ --log ]]; then 
   	   viewlog=`echo ${names[$i]} | sed 's/[-a-zA-Z0-9]*=//'`
       elif [[ "${names[$i]}" =~ -l ]]; then 
           let i=i+1
	   viewlog=${names[$i]}

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

# TO DO:
# rm lock file
# show log files in additional xterm

