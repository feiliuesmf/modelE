#package RegressionTools;
use CommandEntry;
use Env;

my $HYCOM = "E1fzhyc";

my $extraFlags;
$extraFlags{SERIAL}   = "";
$extraFlags{MPI}      = "MPI=YES";
$extraFlags{OPENMP}   = "EXTRA_FFLAGS=-mp MP=YES NPROCS=\$npes";
$extraFlags{SERIALMP} = "EXTRA_FFLAGS=-mp"; # Intel kludge for matching serial and OpenMP.

$extraFlags{E4TcadF40} ="";
$extraFlags{E4arobio_g6cc} ="";
$extraFlags{E4arobio_h4c} ="";
$extraFlags{EM20} ="";
$extraFlags{E4F40} ="";
$extraFlags{E1oM20} ="";
$extraFlags{E_AR5_CADI} ="ESMF=YES";
# Need to set CS variables in the configuration file:
$extraFlags{E4C90L40} ="ESMF=YES FVCUBED=YES FVCUBED_ROOT=/usr/local/other/Fortuna-2_5.noHDF5 MPPDIR=/usr/local/other/Fortuna-2_5.noHDF5/Linux FFTW_ROOT=/discover/nobackup/mkelley5/fftw-3.2.2";

$extraFlags{$HYCOM} ="EXTRA_FFLAGS+=-DCHECK_OCEAN_HYCOM";

$extraFlags{intel} = "";
$extraFlags{gfortran} ="";

# -----------------------------------------------------------------------------
sub createTemporaryCopy 
{
  my $referenceDir = shift;
  my $tempDir = shift;
  my $commandString = "git clone $referenceDir $tempDir";
  print "createTemporaryCopy: $commandString \n";
  return (CommandEntry -> new({COMMAND => $commandString}));
}

# -----------------------------------------------------------------------------
# Build rundeck and compile
sub compileRundeck 
{
  my $env = shift;
  my $installDir = shift;

  my $compiler = $env->{COMPILER};
  my $rundeck  = $env->{RUNDECK};
  my $branch   = $env->{BRANCH};
  my $configuration = $env -> {CONFIGURATION};

  my $resultsDir = $env -> {RESULTS_DIRECTORY};
  $resultsDir .="/$compiler";
 
  my $flags = "$extraFlags{$configuration} $extraFlags{$rundeck} $extraFlags{$compiler}";
  $flags =~ s/(\$npes)/1/eeg;

  my $MODELERC = $env->{MODELERC};

  my $expName;
  if (@_) {$expName = shift}
  else {$expName = "$rundeck.$configuration.$compiler"};

  my $logFile = "$resultsDir/$expName.buildlog";
  unlink($logFile); # delete it

  my $commandString;
  if ($rundeck =~ m/AR5/) 
  { 
    $commandString = <<EOF;
    export MODELERC=$MODELERC;
    cd $installDir/decks;
    make -f new.mk rundeck RUN=$expName RUNSRC=$rundeck;
    make -f new.mk vclean RUN=$expName;
    make -f new.mk -j gcm RUN=$expName $flags COMPILER=$compiler;
EOF
  }
  else 
  {
    $commandString = <<EOF;
    export MODELERC=$MODELERC;
    cd $installDir/decks;
    make rundeck RUN=$expName RUNSRC=$rundeck;
    make vclean RUN=$expName;
    make -j gcm RUN=$expName $flags COMPILER=$compiler;
EOF
  }

  my $binDir = $expName . "_bin";
  print "compileRundeck: $commandString \n";
  return (CommandEntry -> new({COMMAND => $commandString, QUEUE => "", STDOUT_LOG_FILE => "$logFile", COMPILER => $compiler, MODELERC=>$MODELERC, RUNDECK => $rundeck, BRANCH => $branch }));
}

# -----------------------------------------------------------------------------
sub runConfiguration 
{
  my $env = shift;
  my $installDir = shift;
  my $npes = shift;

  my $compiler  = $env->{COMPILER};
  my $rundeck   = $env->{RUNDECK};
  my $branch    = $env->{BRANCH};
  my $configuration = $env->{CONFIGURATION};
  my $resultsDir = $env->{RESULTS_DIRECTORY};
  $resultsDir .="/$compiler";

  print "BRANCH: $branch \n";
  print "resultsDir: $resultsDir \n";

  my $flags = "$extraFlags{$configuration} $extraFlags{$rundeck} $extraFlags{$compiler}";
  $flags =~ s/(\$npes)/$npes/eeg;
  my $expName = "$rundeck.$configuration.$compiler";

  my $suffix;
  my $mpiArgs;
  my $logFile = "$resultsDir/$expName.runlog";

  if ($configuration eq "SERIAL") 
  {
    $suffix = "";
    $mpiArgs = "";
  }
  else 
  {
    $suffix = ".np=$npes";
    $mpiArgs = "-np $npes";
  }
    
  my $MODELERC = $env->{MODELERC};

  my $run1hr;
  my $run1dy;
  if ($rundeck =~ m/AR5/) 
  {
    $run1hr = "export MODELERC=$MODELERC; make -f new.mk rundeck RUN=$expName RUNSRC=$rundeck OVERWRITE=YES; make -f new.mk -j setup RUN=$expName $flags; ../exec/runE_new $expName -np $npes -cold-restart; rm -f $expName/run_status";
    $run1dy = "export MODELERC=$MODELERC; $installDir/exec/editRundeck.sh $expName 48 2 1; make -f new.mk -j setup RUN=$expName $flags; ../exec/runE_new $expName -np $npes -cold-restart; rm -f $expName/run_status";
  }
  else 
  {
    $run1hr = "export MODELERC=$MODELERC; make rundeck RUN=$expName RUNSRC=$rundeck OVERWRITE=YES; make -j setup RUN=$expName $flags; ../exec/runE $expName -np $npes -cold-restart; rm -f $expName/run_status";
    $run1dy = "export MODELERC=$MODELERC; $installDir/exec/editRundeck.sh $expName 48 2 1; make -j setup RUN=$expName $flags; ../exec/runE $expName -np $npes -cold-restart; rm -f $expName/run_status";
  }
  my $restart = "pushd $expName; cp fort.2.nc fort.1.nc ; ./$expName $mpiArgs  ; popd";

  if ($configuration eq "MPI" or $configuration eq "OPENMP") 
  {
    $continue1Day = "./$expName -np $npes ";
  }


  my $commandString = <<EOF;
  export MODELERC=$MODELERC;
  cd $installDir/decks;
  rm -f $expName/fort.2.nc;
  $run1hr;
  if [ -e $expName/fort.2.nc ]; then
    cp $expName/fort.2.nc $resultsDir/$expName.1hr$suffix;
  else
    exit 1;
  fi
  $run1dy;
  if [ -e $expName/fort.1.nc ]; then
    cp $expName/fort.1.nc $resultsDir/$expName.1dy$suffix;
  else
    exit 1;
  fi
  $restart;
  if [ -e $expName/fort.2.nc ]; then
    cp $expName/fort.2.nc $resultsDir/$expName.restart$suffix;
  else
    exit 1;
  fi
EOF
  print "runConfiguration: $commandString \n";
  return (CommandEntry -> new({COMMAND => $commandString, QUEUE => "", STDOUT_LOG_FILE => "$logFile", NUM_PROCS => $npes, COMPILER => $compiler, RUNDECK => $rundeck , BRANCH => $branch }));
}

# -----------------------------------------------------------------------------
sub writeModelErcFile 
{
  my $env = shift;
  my $modelerc = $env->{MODELERC};
  my $commandString .= "mkdir -p $env->{SCRATCH_DIRECTORY}/$env->{COMPILER}\n"; 
  $commandString .= "rm $modelerc\n";
  
  while (my ($var, $value) = each(%$env) ) 
  {
    $commandString .= "echo $var=$value >> $modelerc\n";
  }
  $commandString .= "mkdir -p $env->{DECKS_REPOSITORY} $env->{CMRUNDIR} $env->{SAVEDISK} $env->{EXECDIR} \n";
  print "writeModelErcFile: $commandString \n";
  return (CommandEntry -> new({COMMAND => $commandString}))
}

# -----------------------------------------------------------------------------
sub gitCheckout 
{
  my $env = shift;

  my $scratchDirectory = $env -> {SCRATCH_DIRECTORY};
  my $gitroot          = $env -> {GITROOT};
  my $branch           = $env -> {BRANCH};

  my $commandString = <<EOF;
  pushd $scratchDirectory
  git clone $gitroot $branch
  popd
EOF
  print "gitCheckout: $commandString\n";
  return (CommandEntry -> new({COMMAND => $commandString}))
}

# -----------------------------------------------------------------------------
sub getEnvironment 
{
  my $compiler = shift;
  my $scratchDir = shift;
  my $branch = shift;

  if ($compiler eq "intel") 
  {
    return getIntelEnvironment($scratchDir, $branch);
  }
  else 
  {
    return getGfortranEnvironment($scratchDir, $branch);
  }
}

# -----------------------------------------------------------------------------
sub getIntelEnvironment
{
  my $scratchDir = shift;
  my $branch     = shift;
  print "getIntelEnvironment: $branch \n";

  my $env = {};

  $env->{SCRATCH_DIRECTORY}=$scratchDir;
  $env->{BASELINE_DIRECTORY}=$ENV{MODELEBASELINE}; # NOBACKUP}/modelE_baseline
  $env->{RESULTS_DIRECTORY} =$ENV{REGRESULTS}    ; # NOBACKUP}."/regression_results"
  $env->{DECKS_REPOSITORY}="$scratchDir/decks_repository";
  $env->{CMRUNDIR}="$scratchDir/cmrun";
  $env->{EXECDIR}="$scratchDir/exec";
  $env->{SAVEDISK}="$scratchDir/savedisk";
  $env->{GCMSEARCHPATH}=$ENV{GCMSEARCHPATH}; #discover/nobackup/projects/giss/prod_input_files
  $env->{MP}="no";
  $env->{OVERWRITE}="YES";
  $env->{OUTPUT_TO_FILES}="YES";
  $env->{VERBOSE_OUTPUT}="YES";
  $env->{MPIDISTR}="intel";
  $env->{COMPILER}="intel";
  $env->{GITROOT}=$ENV{MODELROOT}."/$branch";
  if ($branch =~ m/AR5/) 
  {
    $env->{BASELIBDIR}="/usr/local/other_old/esmf/2.2.2rp3_intel-10.1.017_impi-3.2.2.006/Linux";
    $env->{NETCDFHOME}="/usr/local/other/netcdf/3.6.2_intel-10.1.013";
    $env->{PNETCDFHOME}="/discover/nobackup/mkelley5/pnetcdf-1.2.0";
  }
  else 
  {
    $env->{BASELIBDIR5}="/usr/local/other/esmf400rp1/intel11_impi32";
    $env->{NETCDFHOME}="/usr/local/other/netcdf/3.6.2_intel-11.0.083";
    $env->{PNETCDFHOME}="/usr/local/other/pnetcdf/intel11.1.072_impi3.2.2.006";
  }
  $env->{MODELERC}="$scratchDir/intel/modelErc.intel";
  return $env;
}

# -----------------------------------------------------------------------------
sub getGfortranEnvironment 
{
  my $scratchDir = shift;
  my $branch     = shift;
  print "getGFortranvironment: $branch \n";

  my $env = {};
    
  $env->{SCRATCH_DIRECTORY}=$scratchDir;
  $env->{BASELINE_DIRECTORY}=$ENV{MODELEBASELINE}; # NOBACKUP}/modelE_baseline"
  $env->{RESULTS_DIRECTORY}=$ENV{REGRESULTS};      # NOBACKUP}."/regression_results"
  $env->{DECKS_REPOSITORY}="$scratchDir/decks_repository";
  $env->{CMRUNDIR}="$scratchDir/cmrun";
  $env->{EXECDIR}="$scratchDir/exec";
  $env->{SAVEDISK}="$scratchDir/savedisk";
  $env->{GCMSEARCHPATH}=$ENV{GCMSEARCHPATH}; #discover/nobackup/projects/giss/prod_input_files
  $env->{MP}="no";
  $env->{OVERWRITE}="YES";
  $env->{OUTPUT_TO_FILES}="YES";
  $env->{VERBOSE_OUTPUT}="YES";
  $env->{MPIDISTR}="openmpi";
  $env->{COMPILER}="gfortran";
  $env->{GITROOT}=$ENV{MODELROOT}."/$branch"; # NOBACKUP}."/devel/$branch"
  if ($branch =~ m/AR5/) 
  {
    $env->{BASELIBDIR}="/usr/local/other_old/esmf/2.2.2rp3_gcc-4.5_openmpi-1.4.2/Linux";
    $env->{PNETCDFHOME}="/usr/local/other/pnetcdf/gcc4.5_openmpi-1.4.2";
    $env->{PNETCDFHOME}="/discover/nobackup/mkelley5/pnetcdf-1.2.0";
  }
  else 
  {
    $env->{BASELIBDIR5}="/usr/local/other/esmf400rp1/gcc461_openmpi144";
    $env->{NETCDFHOME}="/usr/local/other/netcdf/3.6.2_gcc4.5";
    $env->{PNETCDFHOME}="/usr/local/other/pnetcdf/gcc4.5_openmpi-1.4.2";
  }
  $env->{MODELERC}="$scratchDir/gfortran/modelErc.gfortran";
  return $env;
}

# -----------------------------------------------------------------------------
sub saveForDiffreport()
{
   my $branch = shift;
   my @rundecks;
# Save regTest.cfg in a format that is easily parsed by a bash script
   require 'regTest.cfg';

   if ($branch =~ m/AR5/)
   {
     @rundecks = @AR5decks;
   }
   else
   {
     @rundecks = @decks;
   }
   my $rsize = scalar @rundecks;
   my $csize = scalar @comps;

   my $file = ".regTest.cfg";
   open (FH, "> $file") or die "Can't open $file for write: $!";
   my $i = 0;
   while($i < $rsize)
   {
      print FH "DECK=$rundecks[$i]\n";
      $i++;
   }
   my $i = 0;
   print STDOUT "comps = $csize\n";
   while($i < $csize)
   {
      print FH "COMPILER=$comps[$i]\n";
      $i++;
   }   
   print FH "LEVEL=$level\n";
   print FH "BRANCH=$branch\n";

   close FH or die "Cannot close $file: $!"; 

}

1;
