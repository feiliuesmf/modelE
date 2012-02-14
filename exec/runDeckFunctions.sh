# -------------------------------------------------------------------
# FUNCTIONS
# -------------------------------------------------------------------

# functions used by rundDeck.sh
# History: C.Cruz - created 1/2011 

# -------------------------------------------------------------------
usage() 
# -------------------------------------------------------------------
{
less <<EOF
NAME
     $0 - Script to run modelE decks

DESCRIPTION

     This script runs modelE under various optional settings. By default
     it will run the EM20 rundeck serially and in parallel and perform 
     restart regression. On DISCOVER it will use the Intell 11.1 compiler 
     and Intel MPI. On other machines it will try to use whatever is
     available.

IMPORTANT NOTE

     Code will be retrieved from the specified local git repository 
     (required argument!). Run $0 -h for more information.
    
SYNOPSIS

     $0 <...arguments...>  [..options...]

ARGUMENTS

     -g | --git      : git local repository path
     -h | --help     : this help screen

OPTIONS

     -b | --branch   : git branch (default: master)
     -t | --tmpDir   : location of scratch space (default: $tmpDir/<runid>)
     -m | --modelErc : specified location of .modelErc (default: built-in)
     -d | --deck     : rundeck name (default: EM20)
     -l | --level    : test level (default: debug, 
                       options: gentle, aggressive, insane)
     -s | --serial   : build/run rundeck serially (default: NO) 
     -c | --clean    : clean up scratch space (default: NO)
     -v | --verbose  : verbose (default: YES)
     -r | --restart  : perform restart regression (default: YES)
     -f | --dflags   : use -O0 FFLAGS (default: NO)
     -w | --fvcore   : use CS dycore (default: NO)
     -p | --scmodel  : run single colum model (default: NO, runs serial)

     -a | --endDate  : change the default end date (default: depends on rundeck)
     -o | --endHour  : change the default end hour (default: depends on rundeck)
     -x | --check    : check repository against HEAD on simplex (default: NO)
     -i | --id       : override run ID (default: system generated random number)

     Additionally, on DISCOVER:

     -e | --modEnv   : module environment (default: intel11)
                       where intel11  : comp/intel-11.1.072 
                                        mpi/impi-3.2.2.006
                             intel12  : comp/intel-12.0.1.107 
                                        mpi/impi-3.2.2.006
                             gcc46    : other/comp/gcc-4.6-20110312
					and experimental openmpi 1.4.3
                             gcc461   : experimental gcc4.6.1 
					and experimental openmpi 1.4.3
                             ar5      : comp/intel-10.1.017 
				        mpi/impi-3.2.2.006


EXAMPLE

     runDeck.sh -h 
     runDeck.sh --gitrepo=$HOME/modelE --check=NO
     runDeck.sh -g $HOME/modelE.test --modEnv=gcc46 --modelErc=$HOME/.modelErc.gcc46

     Note that when using short argument names the value of the argument is 
     specified after the argument name separated by a space whereas when using
     using long argument names one must use an equal sign to specify the 
     argument value. Thus the following are equivalent:

     runDeck.sh --gitrepo=$tmpDir/modelE --branch=ar5 --clean=YES --modEnv=gcc45
     runDeck.sh -g $NOBACKUP/modelE -b ar5 -c YES -m gcc45

     One can also mix short and long argument names:

     runDeck.sh --fvcore=YES --deck=E4C90L40 -m=$HOME/.modelErc.cubed -v YES

MORE EXAMPLES

     Run EC12 ( compiled with -O0 ) and change runID to ec12 

     runDeck.sh -g $NOBACKUP/modelE.git -d EC12 -i ec12 -f YES

     Run E4arobio_g6c for 1day cold-restart ( rather than default 1hr setup )

     runDeck.sh -g $HOME/modelE.oceans -d E4arobio_g6c -a 2 -o 0 

     Run E4TcadF40 excluding serial runs and do not compare results with HEAD

     runDeck.sh -g $HOME/modelE.tracers -d E4TcadF40 -s NO -x NO

     Run FVCORE canonical rundeck using specified .modelErc

     runDeck.sh -g $HOME/modelE.fv -w YES -d E4C90L40 -m $HOME/.modelErc.cubed -s NO

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

# Use the smallest possible test level corresposnding to npes=( 1 2 )
   level=debug

# Default is master branch
   gitBranch=

# .modelErc file path - the default file is automatically created here
   rcFile=

# Rundeck name
   rundeck="EM20"
   if [ "$useFVCScore" == "YES" ]; then
      rundeck="E4C90L40"
   fi
# module environment in .modelErc - default is to use Intel 11
   moduleSet="intel11"
   # used in .modelErc:
   COMPILER=intel
   MPIDISTR=intel
   MPIDIR=

# tmpDir is either the nobackup are or $HOME, else use optional argument
   if [[ "$node" =~ borg || "$node" =~ discover || "$node" =~ dali ]]; then
      tmpDir=$NOBACKUP
   else  # not on DISCOVER, e.g. my MAC laptop
      tmpDir=$HOME
   fi

# Run serially - but may want to disable for hi-res rundecks. Anyway all
# runs run in parallel with NPES=1
   runSerial="NO"
   clean="NO"
   verbose="YES"
   restartRegression="YES"

# These are the typical settings for 1-hr setup: DATEE=2,HOURE=1
   endDate=1
   endHour=1
   
# Extra compilation flags
   XFlags="NO"

# Check rundeck changes against those in SIMPLEX repository
   check="NO"

# "diff report" for nc files based on GISS utility. On DISCOVER it is 
# installed and ready to go but user must build it on other systems (and
# make sure it is in the PATH).
   if [[ "$node" =~ borg || "$node" =~ discover || "$node" =~ dali ]]; then
      diffReport=/discover/nobackup/projects/giss/exec/diffreport
   else
      diffReportinPath=`command -v diffreport`
      if [ "$diffReportinPath" == "" ]; then
	 echo " ------ diffreport executable not found..."
	 echo " ------ Will use _diff_ but final report may show a FAILURE."
	 diffReport="diff -q"
      else
         diffReport=$diffReportinPath
      fi
   fi

# fvcore options (default is not to use, else set -fvcore option)
   FVCORE=
   FVCUBED=
   FVCORE_ROOT=
   FVCUBED_ROOT=
   MPPDIR=
   FFTW_ROOT=

}

# -------------------------------------------------------------------
initEnvironment()
# -------------------------------------------------------------------
{

# tmpDir is either $NOBACKUP (DISCOVER), $HOME (not DISCOVER) or custom
# Use random run id, unless specified in command line 
   if [ "$runID" == "" ]; then
      runID=$RANDOM.$$
      workPath=$tmpDir/tmp.$runID
   else
      workPath=$tmpDir/$runID
   fi

# Checkpoint initial path
   startDirectory=`pwd`

# scratch space
   ( umask 022 && ( mkdir -p $workPath > /dev/null 2>&1 ) ) || {
        echo " ~~~ Could not create temporary directory $workPath" 1>&2
        exit 1
   }
   cd $workPath

# set email address and log file names
   setLogNames

   diagMessage " *** Setup $rundeck testing environment ***"

# check/set module environment
   checkSetModEnv "$moduleSet"

# get code from repository
   gitClone

# create a default .modelErc file that works with development branch
   createRcFile

# set test level: gentle (3 lats per proc), aggressive (2 lats), insane
   setTestLevel

# Set EXTRA_FFLAGS 
   if [ "$XFlags" == "YES" ]; then
      EXTRA_FFLAGS="-O0"
   fi

# exclude PFUNIT
   unset PFUNIT

   if [ "$verbose" == "YES" ]; then
      printInfo
   fi
}

# -------------------------------------------------------------------
setLogNames()
# -------------------------------------------------------------------
{
# log files
   toEmail=${workPath}/"email-modele-log"
   makeLog=${workPath}/"make-modele-log"
   rm -f $toEmail $makeLog

# On Linux (DISCOVER) get email from environment
   if [[ "$node" =~ borg ]]; then
#       emailAddress=$USER@nccs.nasa.gov
       emailAddress=ccruz@nccs.nasa.gov
   else
       emailAddress=$USER@localhost
   fi
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
         diagMessage " ------ MPI compilation will FAIL"
      else
         diagMessage " ------ MPIRUN : "$mpiCheck
         mpiCheck=`mpif90 --version | head -1`
         diagMessage " ------ MPIF90 will use the following compiler : "$mpiCheck
      fi  
      gitCheck=`which git`
      if [ "$gitCheck" == "" ]; then
         finalize 1 " --- No GIT installation found."
      fi
   fi
}

# -------------------------------------------------------------------
gitClone()
# -------------------------------------------------------------------
{
# Two git repositories are required. The git repository required from
# the command line is henceforth referred simply as the gitRepository.
# The changes we make in that repository are compared against the baseline
# repository, henceforth referred as the gitBaseRepository, typically 
# in the simplex server. Since we are unable to ssh from compute nodes then 
# our gitBaseRepository must be "local" and, unfortunately, hardwired to 
# this script. When we find a bway to ssh from compute nodes then:
#
#    gitBaseRepository=simplex.giss.nasa.gov:/giss/gitrepo/modelE.git
#
   diagMessage " --- Cloning repositories..."
   if [ "$check" == "YES" ]; then
      if [[ "$node" =~ borg ]]; then
         gitBaseRepository=/discover/nobackup/ccruz/devel/modelE.clones/simplex
      else
         gitBaseRepository=$HOME/models/devel/modelE.clones
      fi
      # first get the gitBaseRepository:
      if [ "$gitBranch" == "" ]; then
         git clone --quiet $gitBaseRepository modelE.baseline > /dev/null 2>&1
      else
         git clone --quiet -b $gitBranch $gitBaseRepository modelE.baseline > /dev/null 2>&1
      fi

      if [ "$?" -ne "0" ]; then
         clean="YES"
         finalize 1 " --- Failed to clone git repository: $gitBaseRepository"
      fi
   fi

   # then get the git experimental repository:
   if [ "$gitBranch" == "" ]; then
      git clone --quiet $gitRepository modelE > /dev/null 2>&1
   else
      git clone --quiet -b $gitBranch $gitRepository modelE > /dev/null 2>&1
   fi
   if [ "$?" -ne "0" ]; then
      clean="YES"
      finalize 1 " --- Failed to clone git repository: $gitRepository"
   fi
   diagMessage " --- Done cloning."
}

# -------------------------------------------------------------------
setModules()
# -------------------------------------------------------------------
{
# Set module environment based  on input argument. Here we also set
# the library paths used in .modelErc

    local moduleSet=$1

    source /usr/share/modules/init/bash
    module purge
    # NOTE: DEFAULT is intel11
    if [ "$moduleSet" == "intel11" ]; then
       module load comp/intel-11.1.072 mpi/impi-3.2.2.006
       PNETCDFHOME=/discover/nobackup/mkelley5/pnetcdf-1.2.0
       NETCDFHOME=/discover/nobackup/projects/gmao/share/dao_ops/Baselibs/v3.2.0_build3/Linux
       ESMFHOME=/discover/nobackup/projects/gmao/share/dao_ops/Baselibs/v3.2.0_build3/Linux
    elif [ "$moduleSet" == "intel12" ]; then
       module load comp/intel-12.0.1.107 mpi/impi-3.2.2.006
       PNETCDFHOME=/usr/local/other/pnetcdf/intel12.0.1.107_impi3.2.2.006
       NETCDFHOME=/usr/local/other/netcdf/3.6.2_intel-12.0.1.107
       ESMFHOME=/usr/local/other/esmf5/intel12.0.1.107_impi3.2.2.006/Linux
    elif [ "$moduleSet" == "gcc46" ]; then
       module load other/comp/gcc-4.6-20110312
       PNETCDFHOME=/usr/local/other/pnetcdf/gcc4.5_openmpi-1.4.2
       NETCDFHOME=/usr/local/other/netcdf/3.6.2_gcc4.6
       ESMFHOME=/usr/local/other/esmf5/gcc4.5_openmpi-1.4.2/Linux
       COMPILER=gfortran
       MPIDISTR=openmpi
       MPIDIR=/gpfsm/dnb32/ccruz/Baselibs/openmpi/1.4.3-gcc-4.6
       # need these two exports until this openmpi is a module
       export LD_LIBRARY_PATH=/gpfsm/dnb32/ccruz/Baselibs/openmpi/1.4.3-gcc-4.6/lib:${LD_LIBRARY_PATH}
       export PATH=/gpfsm/dnb32/ccruz/Baselibs/openmpi/1.4.3-gcc-4.6/bin:${PATH}
    elif [ "$moduleSet" == "gcc461" ]; then
       module load other/comp/gcc-4.6.1-RC-20110620
       PNETCDFHOME=/usr/local/other/pnetcdf/gcc4.5_openmpi-1.4.2
       NETCDFHOME=/usr/local/other/netcdf/3.6.2_gcc4.6
       ESMFHOME=/usr/local/other/esmf5/gcc4.5_openmpi-1.4.2/Linux
       COMPILER=gfortran
       MPIDISTR=openmpi
       MPIDIR=/gpfsm/dnb32/ccruz/Baselibs/openmpi/1.4.3-gcc-4.6.1
       export LD_LIBRARY_PATH=/gpfsm/dnb32/ccruz/Baselibs/openmpi/1.4.3-gcc-4.6.1/lib:${LD_LIBRARY_PATH}
       export PATH=/gpfsm/dnb32/ccruz/Baselibs/openmpi/1.4.3-gcc-4.6.1/bin:${PATH}
    elif [ "$moduleSet" == "ar5" ]; then
       module load comp/intel-10.1.017 mpi/impi-3.2.2.006
       PNETCDFHOME=/discover/nobackup/mkelley5/pnetcdf-1.2.0
       NETCDFHOME=/usr/local/other/netcdf/3.6.2_intel-10.1.013
       ESMFHOME=/usr/local/other/baselibs/ESMF222rp3_NetCDF362b6_10.1.017_intelmpi/Linux
    else
      finalize 1 " --- Module set $moduleSet is not available."
    fi
    # finally load GIT module
    module load other/git-1.7.3.4
}

# -------------------------------------------------------------------
createRcFile()
# -------------------------------------------------------------------
{
# If modelErc is not specified we will use the following modelErc:

   if [ "$rcFile" == "" ]; then

# This rc file should work with the development branch and assumes 
# FVCORE=NO

# FVCS core
   if [ "$useFVCScore" == "YES" ]; then
      if [ "$moduleSet" != "intel11" ]; then
          diagMessage " --- FVCS core was not build with this $moduleSet."
          diagMessage " --- *** runDeck.sh will probably fail *** "
      fi
      FVCORE=YES
      FVCUBED=YES
      # temporary
      FVCUBED_ROOT=/usr/local/other/Fortuna-2_5.noHDF5
      MPPDIR=/usr/local/other/Fortuna-2_5.noHDF5/Linux
      FFTW_ROOT=/discover/nobackup/mkelley5/fftw-3.2.2
   fi

cat << EOF > .modelErc
# directories and general option
DECKS_REPOSITORY=$workPath/decks
CMRUNDIR=$workPath/cmrun
EXECDIR=$workPath/exec
SAVEDISK=$workPath/out
GCMSEARCHPATH=/discover/nobackup/projects/giss/prod_input_files
OUTPUT_TO_FILES=YES
VERBOSE_OUTPUT=YES

# netcdf
NETCDFHOME=$NETCDFHOME
PNETCDFHOME=$PNETCDFHOME

# compiler
COMPILER=$COMPILER

# mpi
MPIDISTR=$MPIDISTR
MPIDIR=$MPIDIR

MP=NO

EOF

if [ "$moduleSet" == "ar5" ]; then
cat << EOF >> .modelErc
BASELIBDIR=$ESMFHOME
EOF
else
cat << EOF >> .modelErc
BASELIBDIR5=$ESMFHOME
EOF
fi

cat << EOF >> .modelErc

# fvcore options
FVCORE=$FVCORE
FVCUBED=$FVCUBED
FVCUBED_ROOT=$FVCUBED_ROOT
MPPDIR=$MPPDIR

# other options
FFTW_ROOT=$FFTW_ROOT

EOF

      rcFile=$workPath/.modelErc

   else

      if [ ! -e "$rcFile" ]; then
         finalize 1 " --- MODELERC=$rcFile does not exist."
      fi
      
# Since we are using a specific .modelErc we copy it to $workPath
# and edit SAVEDISK and DECKS_REPOSITORY to exists under $workpath
# EXECDIR and CMRUNDIR are not being used so ther are not modified.

      cp $rcFile $workPath/.modelErc
      out="SAVEDISK=$workPath/out"
      decks="DECKS_REPOSITORY=$workPath/decks"
      cmrun="CMRUNDIR=$workPath/cmrun"
      exec="EXECDIR=$workPath/exec"
      while read line
      do
         if [[ "$line" =~ SAVEDISK= ]]; then
            echo  $out >> $workPath/.foo  
         elif [[ "$line" =~ DECKS_REPOSITORY= ]]; then
            echo  $decks >> $workPath/.foo  
         elif [[ "$line" =~ CMRUNDIR= ]]; then
            echo  $cmrun >> $workPath/.foo  
         elif [[ "$line" =~ EXECDIR= ]]; then
            echo  $exec >> $workPath/.foo  
         else
            echo $line >> $workPath/.foo  
         fi
      done < $workPath/.modelErc
      mv $workPath/.foo $workPath/.modelErc
      rcFile=$workPath/.modelErc

   fi

   if [ ! -d $workPath/decks ]; then mkdir $workPath/decks; fi
   if [ ! -d $workPath/out   ]; then mkdir $workPath/out;   fi
   if [ ! -d $workPath/cmrun ]; then mkdir $workPath/cmrun; fi
   if [ ! -d $workPath/exec  ]; then mkdir $workPath/exec;  fi

}

# -------------------------------------------------------------------
setTestLevel()
# -------------------------------------------------------------------
{
# In this function we determine, based on the "level" setting, how many
# cpus to use when running a particular rundeck. 

   if [[ "$level" =~ gentle ]]; then
      # EM20 = 1,15
      # E4F40, E4TcadF40, E4arobio_H4c, E4arobio_g6c = 1,4,30
      if [[ "$rundeck" =~ EM20 ]]; then
         npes=( 1 15 )
      elif [[ "$rundeck" =~ E4Tcad || "$rundeck" =~ E4arobio  ]]; then
         npes=( 23 )
      elif [[ "$rundeck" =~ E4C90 ]]; then
         npes=( 6 )
      else 
         npes=( 1 6 )
      fi
   elif [[ "$level" =~ aggressive ]]; then
      # EM20 = 1,4,23
      # E4F40 = 1,4,45
      # E4TcadF40, E4arobio_H4c, E4arobio_g6c = 1,4,30
      if [[ "$rundeck" =~ EM20 ]]; then
         npes=( 1 23 )
      elif [[ "$rundeck" =~ E4F40 ]]; then
         npes=( 1 45 )
      elif [[ "$rundeck" =~ E4Tcad || "$rundeck" =~ E4arobio ]]; then
         npes=( 45 )
      elif [[ "$rundeck" =~ E4C90 ]]; then
         npes=( 6 48 )
      else 
         npes=( 1 48)
      fi
   elif [[ "$level" =~ insane ]]; then
      # EM20 = 1,44
      # E4F40, E4TcadF40, E4arobio_g6c = 1,88
      # E4arobio_H4c = 1,44
      if [[ "$rundeck" =~ EM20 ]]; then
         npes=( 1 44 )
      elif [[ "$rundeck" =~ E4F40 || "$rundeck" =~ E4Tcad || "$rundeck" =~ E4arobio_g6c ]]; then
         npes=( 88 )
      elif [[ "$rundeck" =~ E4arobio_h4c ]]; then
         npes=( 44 )
      elif [[ "$rundeck" =~ E4C90 ]]; then
         npes=( 6 84 )
      else 
         npes=( 1 84)
      fi
   else
      npes=( 1 2 )
   fi
}

# -------------------------------------------------------------------
buildAndRun() 
# -------------------------------------------------------------------
{
   local val
   local cnt
   local i
   local codebase

   diagMessage " *** Running $rundeck ... please wait ***"
   diagMessage " ------------------------------------------"
   if [ "$check" == "YES" ]; then
      codebase=base
      diagMessage " --- BASELINE ---"
      cd $workPath/modelE.baseline/decks
      buildAndRunCodeBase "$codebase"
   fi
   codebase=exp
   diagMessage " --- EXPERIMENTAL ---"
   cd $workPath/modelE/decks
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

   if [[ "$runSerial" == "YES" || "$check" == "YES" || "$runSCM" == "YES" ]]; then
      buildSerial "${rundeck}.serial.${codebase}"
      runModel "${rundeck}.serial.${codebase}" 0

# clean up for MPI build
      diagMessage " ------ gmake vclean"
      gmake vclean 1>> $makeLog 2>&1
      wait
   fi

   if [[ "$codebase" != "base" || "$runSCM" == "NO" ]]; then
      cnt=0
      buildOption=YES
      for val in "${npes[@]}"; do
	 if [ $cnt -gt 0 ]; then
            buildOption="NO"
	 fi 
	 buildMpi "${rundeck}.${codebase}" $val "$buildOption"
	 let cnt=$cnt+1
# Now run modelE setup in parallel with NPES=val
	 if [[ $endDate -ne 1 || $endHour -ne 1 ]]; then
            editRundeck "$rundeck.${codebase}.R" 960 $endDate $endHour
	 fi
	 runModel "${rundeck}.${codebase}" "$val"
	 saveState "${rundeck}.${codebase}" "$val" "NO"
# Run restart regression (compares 25 hr restart vs 24+1hr restart)
	 if [ "$restartRegression" == "YES" ]; then
            restartRegression "${rundeck}.${codebase}" "$val"
	 fi
      done
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

# make gcm
   echo " ====== START SERIAL BUILD ====== " >> $makeLog
   diagMessage " ------ gmake -j gcm RUN=$deck"
   if [ "$EXTRA_FFLAGS" == "" ]; then   
      gmake -j gcm RUN=$deck MODELERC=$rcFile 1>> $makeLog 2>&1
   else
      gmake -j gcm RUN=$deck MODELERC=$rcFile EXTRA_FFLAGS=$EXTRA_FFLAGS 1>> $makeLog 2>&1
   fi
   wait
   exitIfNofile "$deck.exe"
}

# -------------------------------------------------------------------
buildMpi() 
# -------------------------------------------------------------------
{
   local deck=$1
   local npe=$2
   local buildOpt=$3

   if [ "$buildOpt" == "YES" ]; then

      diagMessage " --- Building modelE MPI version"

# make rundeck
      makeRundeck "$deck"
      exitIfNofile "$deck.R"

# make gcm
      echo " ====== START MPI BUILD ====== " >> $makeLog
      if [ "$useFVCScore" == "YES" ]; then # we need ESMF
         diagMessage " ------ gmake -j gcm RUN=$deck ESMF=YES" 
         if [ "$EXTRA_FFLAGS" == "" ]; then   
            gmake -j gcm RUN=$deck MODELERC=$rcFile ESMF=YES 1>> $makeLog 2>&1
         else
            gmake -j gcm RUN=$deck MODELERC=$rcFile ESMF=YES EXTRA_FFLAGS=$EXTRA_FFLAGS 1>> $makeLog 2>&1
         fi
      else # we use MPI
         diagMessage " ------ gmake -j gcm RUN=$deck MPI=YES" 
         if [ "$EXTRA_FFLAGS" == "" ]; then   
            gmake -j gcm RUN=$deck MODELERC=$rcFile MPI=YES 1>> $makeLog 2>&1
         else
            gmake -j gcm RUN=$deck MODELERC=$rcFile MPI=YES EXTRA_FFLAGS=$EXTRA_FFLAGS 1>> $makeLog 2>&1
         fi
      fi
      wait
      exitIfNofile "$deck.exe"

   fi # end if buildOpt
}

# -------------------------------------------------------------------
makeRundeck() 
# -------------------------------------------------------------------
{ 
   local deck=$1
   if [ ! -e ../templates/$rundeck.R ]; then
      finalize 1 " --- RUNSRC=$rundeck does not exist."
   fi
   diagMessage " ------ gmake rundeck RUN=$deck RUNSRC=$rundeck"
   gmake rundeck RUN=$deck RUNSRC=$rundeck OVERWRITE=YES MODELERC=$rcFile 1>> $makeLog 2>&1
}

# -------------------------------------------------------------------
runModel() 
# -------------------------------------------------------------------
{   
   local deck=$1
   local ncpu=$2
   if [ $ncpu -eq 0 ]; then
      diagMessage " --- Run $deck (serial mode)"
      gmake -j setup RUN=$deck MODELERC=$rcFile 1>> $makeLog 2>&1
      MODELERC=$rcFile ../exec/runE $deck -np 1 -cold-restart 1>> $makeLog 2>&1
   else
       diagMessage " --- Run $deck (MPI mode,  NPES=$ncpu)"
       gmake -j setup RUN=$deck MODELERC=$rcFile MPI=YES 1>> $makeLog 2>&1
       MODELERC=$rcFile ../exec/runE $deck -np $ncpu -cold-restart 1>> $makeLog 2>&1
   fi
   checkSetup "$deck"
}

# -------------------------------------------------------------------
checkSetup() 
# -------------------------------------------------------------------
{   
   local deck=$1
   if [ -e $workPath/out/$deck/run_status ]; then
      rcCode=`head -1 $workPath/out/$deck/run_status`
      rcMess=`tail -1 $workPath/out/$deck/run_status`
      diagMessage " ------ run_status : rc=$rcCode; $rcMess"
      if [ $rcCode -ne 13 ] && [ $rcCode -ne 12 ] ; then
         finalize 1 " ------ Problem encountered while running $deck";
      fi
   else
      finalize 1 " ------ $workPath/out/$deck/run_status does not exist";
   fi 
}

# -------------------------------------------------------------------
compareRuns() 
# -------------------------------------------------------------------
{
   local n
   local cnt
   local rc
   local out01=$rundeck.serial.base
   local out02=$rundeck.serial.exp
   local out2=$rundeck.exp
   declare -a id
   id=( )

   cd $workPath/out

   diagMessage " ------------------------------------------"
   cnt=0
   for n in "${npes[@]}"; do
      id[${#id[*]}]="NPES${n}"
      let cnt=$cnt+1
   done

   diagMessage " --- DIFF results: Strong reproducibility"   
   if [ "$runSerial" == "YES" ]; then
      for (( n=0 ;  n<${#npes[*]}; n++ )); do    
         rc=`$diffReport $out02/fort.2.nc $out2/fort.2.nc.${id[$n]} > .sz`
         sz=`ls -la .sz | awk '{print $5}'`
         checkStatus $sz " ------ SERIAL vs ${id[$n]}: "
         rm -f .sz
      done
   fi

   let cnt=${#npes[*]}-1
   for (( n=0 ;  n<$cnt; n++ )); do
      let n1=$n+1
      rc=`$diffReport $out2/fort.2.nc.${id[0]} $out2/fort.2.nc.${id[$n1]} > .sz`
      sz=`ls -la .sz | awk '{print $5}'`
      checkStatus $sz " ------ ${id[0]} vs ${id[$n1]}: "
      rm -f .sz
   done

   if [ "$check" == "YES" ]; then
      diagMessage " --- DIFF results: Has model changed?"   
      if [ "$runSerial" == "YES" ]; then
         rc=`$diffReport $out01/fort.2.nc $out02/fort.2.nc > .sz`
         sz=`ls -la .sz | awk '{print $5}'`
         checkStatus $sz " ------ BASELINE (SERIAL) vs EXP (SERIAL): "
         rm -f .sz
      fi

      let cnt=${#npes[*]}
      for (( n=0 ;  n<$cnt; n++ )); do
         rc=`$diffReport $out01/fort.2.nc $out2/fort.2.nc.${id[$n]} > .sz`
         sz=`ls -la .sz | awk '{print $5}'`
         checkStatus $sz " ------ BASELINE (SERIAL) vs EXP (${id[$n]}): "
         rm -f .sz
      done
   fi
   unset id
   diagMessage " ------------------------------------------"
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

   cd $workPath/out

   cnt=0
   for n in "${npes[@]}"; do
      id[${#id[*]}]="NPES${n}"
      let cnt=$cnt+1
   done

   diagMessage " --- DIFF results - Restart regression"   

   let cnt=${#npes[*]}
   for (( n=0 ;  n<$cnt; n++ )); do
      rc=`$diffReport $out2/restart.1dy.${id[$n]} $out2/fort.2.nc > .sz`
      sz=`ls -la .sz | awk '{print $5}'`
      checkStatus $sz " ------ restart regression for ${id[$n]} : "
      rm -f .sz
   done
   unset id
   diagMessage " ------------------------------------------"
}

# -------------------------------------------------------------------
restartRegression() 
# -------------------------------------------------------------------
{ 
   local deck=$1
   local ncpu=$2
   local regDeck=${deck}.reg
   diagMessage " --- Regression runs"
# Make a new rundeck - it will have the .reg suffix
   makeRundeck "$regDeck"
# Edit rundeck to run 1 hr past day 2 (any assumptions here?) and
# create a checkpoint (fort.1) at day 2.
# Note that 48 = # time steps in one day
   editRundeck "$regDeck.R" 48 2 1
# Run model from t0 -> t2 (running one day)
   runModel "$regDeck" "$ncpu"   
# Remove run_status - this is only a problem on discover (why?)
   rm -f $workPath/out/$regDeck/run_status
# Save t2 (fort.2) into t1 (fort.1) , and keep a copy of t2 for
# comparison with t2'
   saveState "$regDeck" "$ncpu" "YES"
   cd $workPath/out/$regDeck
# Run model from from t1 -> t2'   
   diagMessage " --------- ./$regDeck -np $ncpu" 
   ./$regDeck -np "$ncpu" 1>> $makeLog 2>&1     
   cd $workPath/modelE/decks
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

# this section may be omitted, once all old style rundecks are modified
   a=$( grep ^' '${eof2} $deck )
   if [[ $a = '' ]];then
      eof2='&END'
   fi

   a=$( grep -n ${eof1}     $deck | head -1 ) ; n1=${a%%:*}
   a=$( grep -n ^' '${eof2} $deck | head -1 ) ; n2=${a%%:*}

   cp ${deck} templ
   head -$(( n1-1 )) templ                   > ${deck}
   echo "${ndisk_line}"                      >> ${deck}
   tail +${n1} templ | head -$(( n2 - n1 ))  >> ${deck}
   echo "${end_hour_line}"                   >> ${deck}
   echo " ISTART=2, ${end_hour_line}"         >> ${deck}
   tail +${n2} templ                         >> ${deck}
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
      cp $workPath/out/$deck/fort.1.nc $workPath/out/$deck/restart.1dy.$id
      cp $workPath/out/$deck/fort.2.nc $workPath/out/$deck/fort.1.nc
   else
      cp $workPath/out/$deck/fort.2.nc $workPath/out/$deck/fort.2.nc.$id
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
      message " @@@ $rundeck test: FAILURE @@@"   
      message "$exitMessage"
   fi
   diagMessage " *** $rundeck testing is complete. ***"
# Only e-mail results in PBS_BATCH mode.
   subject="$0 $rundeck $moduleSet results"
   if [[ "$pbsType" =~ BATCH ]]; then
      if [ $exitCode -eq 0 ]; then
         mail -s " $rundeck $moduleSet test: SUCCESS" $emailAddress
      else
         mail -s "$subject" $emailAddress < $toEmail
      fi
# Write results to a file and then email from within the calling script
#      TEMPNAME="$PBS_O_WORKDIR/$level.results"
#      diagMessage "write to $TEMPNAME"
#      cat $toEmail >> $TEMPNAME
   fi
   clean
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
      diagMessage "$msg""0-diff"
      return
   fi
   message "$msg"" NON-ZERO diff <------ FAILURE"
}

# -------------------------------------------------------------------
exitIfNofile() 
# -------------------------------------------------------------------
{ 
   local file=$1
   result=`find . -name $file`
   if [ "$result" == "" ]; then
      finalize 1 " ------ $file does not exist. Please see $makeLog"
   fi
}

# -------------------------------------------------------------------
warnIfNofile() 
# -------------------------------------------------------------------
{ 
   local file=$1
   result=`find . -name $file`
   if [ "$result" == "" ]; then
      diagMessage " ------ $file does not exist."
   fi
}

# -------------------------------------------------------------------
clean() 
# -------------------------------------------------------------------
{
   diagMessage " --- Clean up."
   if [ "$clean" == "YES" ]; then
      rmWorkPath
   else
      diagMessage " ------ Not removing $workPath" 
   fi
   diagMessage " --- Done."
}

# -------------------------------------------------------------------
rmWorkPath() 
# -------------------------------------------------------------------
{
   echo " ------ Removing $workPath" 
   # Just in case user decided to change the code...
   if [[ "$workPath" == "$HOME" || "$workPath" == "$NOBACKUP" ]]; then
      echo " ~~~ CLEAN aborted: cannot remove $workPath"; exit 1
   fi
   rm -rf $workPath
}

# -------------------------------------------------------------------
printInfo()
# -------------------------------------------------------------------
{   
   diagMessage " *** Options ***"
   diagMessage "    EXPERIMENTAL GIT repository: $gitRepository"
   if [ "$check" == "YES" ]; then
      diagMessage "    BASELINE GIT repository: $gitBaseRepository"
   fi
   if [ "$gitBranch" != "" ]; then
      diagMessage "      (Will retrieve from git branch=$gitBranch)"
   fi
   if [[ "$pbsType" =~ BATCH ]]; then
      diagMessage "    moduleEnv  = $moduleSet"
   fi
   diagMessage "    runID      = $runID"
   diagMessage "    tmpDir     = $workPath"
   diagMessage "    rcFile     = $rcFile"
   diagMessage "    runDeck    = $rundeck"
   diagMessage "    runSerial  = $runSerial <- run serial cases"
   diagMessage "    restartReg = $restartRegression <- restart regression"
   diagMessage "    checkBase  = $check <- check against base repository"
   diagMessage "    cleanTmp   = $clean <- clean up scratch space when done"
   diagMessage "    extraFlags = $XFlags <- Use -O0 compilation"
   diagMessage "    testLevel  = $level"
   diagMessage "    endDate    = $endDate"
   diagMessage "    endHour    = $endHour"
}

# -------------------------------------------------------------------
diagMessage()
# -------------------------------------------------------------------
{
   local message=$1
   if [ "$verbose" == "YES" ]; then
      echo "$message"
   fi
}

# -------------------------------------------------------------------
message()
# -------------------------------------------------------------------
{
   local message=$1
   echo "$message"
   echo "$message" >> $toEmail
}

# -------------------------------------------------------------------
parseOptions()
# -------------------------------------------------------------------
{
# Edit this function only if adding a new option.

# first argument MUST be -g or --git or -h or --help
# Eventually -g (or --git) may not be necessary as we will clone directly
# from the simplex repository.

   i1=0
   if [[ "${names[i1]}" =~ --git ]]; then 
      gitRepository=`echo ${names[i1]} | sed 's/[-a-zA-Z0-9]*=//'`
   elif [[ "${names[i1]}" == "-g" ]]; then 
      let i1=i1+1
      gitRepository=${names[i1]}
      if [[ "$gitRepository" == "" || "${gitRepository[@]:0:1}" == "-" ]]; then
         echo " ### Blank name for git repository"; exit 1
      fi
   elif [[ "${names[0]}" =~ --help ]]; then 
      usage; exit 0
   elif [[ "${names[0]}" =~ -h ]]; then 
      usage; exit 0
   else
      echo " ### Unknown first argument ${names[0]}"; exit 1
   fi
   let i1=i1+1

# Options...lots of them
   for (( i = i1 ; i < ${#names[@]} ; i++ ))
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

      elif [[ "${names[$i]}" =~ --tmpDir ]]; then 
	 tmpDir=`echo ${names[$i]} | sed 's/[-a-zA-Z0-9]*=//'`
      elif [[ "${names[$i]}" =~ -t ]]; then 
         let i=i+1
	 tmpDir=${names[$i]}

      elif [[ "${names[$i]}" =~ --branch ]]; then 
	 gitBranch=`echo ${names[$i]} | sed 's/[-a-zA-Z0-9]*=//'`
      elif [[ "${names[$i]}" =~ -b ]]; then 
         let i=i+1
	 gitBranch=${names[$i]}

      elif [[ "${names[$i]}" =~ --modelErc ]]; then 
	 rcFile=`echo ${names[$i]} | sed 's/[-a-zA-Z0-9]*=//'`
      elif [[ "${names[$i]}" =~ -m ]]; then 
         let i=i+1
	 rcFile=${names[$i]}

      elif [[ "${names[$i]}" =~ --serial ]]; then 
	 runSerial=`echo ${names[$i]} | sed 's/[-a-zA-Z0-9]*=//'`
      elif [[ "${names[$i]}" =~ -s ]]; then 
         let i=i+1
	 runSerial=${names[$i]}

      elif [[ "${names[$i]}" =~ --clean ]]; then 
	 clean=`echo ${names[$i]} | sed 's/[-a-zA-Z0-9]*=//'`
      elif [[ "${names[$i]}" =~ -c ]]; then 
         let i=i+1
         clean=${names[$i]}

      elif [[ "${names[$i]}" =~ --restart ]]; then 
	 restartRegression=`echo ${names[$i]} | sed 's/[-a-zA-Z0-9]*=//'`
      elif [[ "${names[$i]}" =~ -r ]]; then 
         let i=i+1
	 restartRegression=${names[$i]}

      elif [[ "${names[$i]}" =~ --dflags ]]; then 
	 XFlags=`echo ${names[$i]} | sed 's/[-a-zA-Z0-9]*=//'`
      elif [[ "${names[$i]}" =~ -f ]]; then 
         let i=i+1
	 XFlags=${names[$i]}

      elif [[ "${names[$i]}" =~ --fvcore ]]; then 
	 useFVCScore=`echo ${names[$i]} | sed 's/[-a-zA-Z0-9]*=//'`
      elif [[ "${names[$i]}" =~ -w ]]; then 
         let i=i+1
	 useFVCScore=${names[$i]}

      elif [[ "${names[$i]}" =~ --scmodel ]]; then 
	 runSCM=`echo ${names[$i]} | sed 's/[-a-zA-Z0-9]*=//'`
      elif [[ "${names[$i]}" =~ -p ]]; then 
         let i=i+1
	 runSCM=${names[$i]}

      elif [[ "${names[$i]}" =~ --verbose ]]; then 
	 verbose=`echo ${names[$i]} | sed 's/[-a-zA-Z0-9]*=//'`
      elif [[ "${names[$i]}" =~ -v ]]; then 
         let i=i+1
	 verbose=${names[$i]}

      elif [[ "${names[$i]}" =~ --endDate ]]; then 
	 endDate=`echo ${names[$i]} | sed 's/[-a-zA-Z0-9]*=//'`
      elif [[ "${names[$i]}" =~ -a ]]; then 
         let i=i+1
	 endDate=${names[$i]}

      elif [[ "${names[$i]}" =~ --endHour ]]; then 
	 endHour=`echo ${names[$i]} | sed 's/[-a-zA-Z0-9]*=//'`
      elif [[ "${names[$i]}" =~ -o ]]; then 
         let i=i+1
	 endHour=${names[$i]}

      elif [[ "${names[$i]}" =~ --level ]]; then 
	 level=`echo ${names[$i]} | sed 's/[-a-zA-Z0-9]*=//'`
      elif [[ "${names[$i]}" =~ -l ]]; then 
         let i=i+1
	 level=${names[$i]}

      elif [[ "${names[$i]}" =~ --id ]]; then 
	 runID=`echo ${names[$i]} | sed 's/[-a-zA-Z0-9]*=//'`
      elif [[ "${names[$i]}" =~ -i ]]; then 
         let i=i+1
	 runID=${names[$i]}

      elif [[ "${names[$i]}" =~ --check ]]; then 
	 check=`echo ${names[$i]} | sed 's/[-a-zA-Z0-9]*=//'`
      elif [[ "${names[$i]}" =~ -x ]]; then 
         let i=i+1
	 check=${names[$i]}

      elif [[ "${names[$i]}" =~ --help ]]; then 
         usage; exit 0
      elif [[ "${names[$i]}" =~ -h ]]; then 
         usage; exit 0

      else
         echo " ### Unknown argument ${names[$i]}"; exit 1
      fi
   done
}

# Do not edit the following function...unless you intend to replace 
# or improve the command line parser.
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
    rmWorkPath
    cleanexit $exit_status
} 

# -------------------------------------------------------------------
tableStatus() 
# -------------------------------------------------------------------
{
# Under construction
   message " --- Test results"   
   line="   0 1 3 6 8  "
   line1=" +-----------+"
   output[1]="0| .         |"
   output[2]="1|   .       |"
   output[3]="3|     .     |"
   output[4]="6|       .   |"
   output[5]="8|         . |"

#tput bold   # Bold print.
   echo -e "\033[1m"    # Bold. May not be portable?
   echo -e "\t$line"
   echo -e "\t$line1"
   echo -e "\t${output[1]}"
   echo -e "\t${output[2]}"
   echo -e "\t${output[3]}"
   echo -e "\t${output[4]}"
   echo -e "\t${output[5]}"
   echo -e "\t$line1"
   echo -e "\n"
#tput sgr0   # Reset terminal.
   echo -e "\033[0m"    # Turn off bold.
}

