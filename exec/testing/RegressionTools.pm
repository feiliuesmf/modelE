use CommandEntry;
use Env;

my $extraFlags;
$extraFlags{SERIAL}   = "";
$extraFlags{MPI}      = "MPI=YES";
$extraFlags{OPENMP}   = "EXTRA_FFLAGS=-mp MP=YES NPROCS=\$npes";
$extraFlags{SERIALMP} = "EXTRA_FFLAGS=-mp"; # Intel kludge for matching serial and OpenMP.

$extraFlags{E4TcadF40} ="";
$extraFlags{E4arobio_g6c} ="";
$extraFlags{E4arobio_h4c} ="";
$extraFlags{EM20} ="";
$extraFlags{E4F40} ="";
$extraFlags{E1oM20} ="";
$extraFlags{E_AR5_CADI} ="ESMF=YES";
# Need to set CS variables in the configuration file:
$extraFlags{E4C90L40} ="ESMF=YES FVCUBED=YES FVCUBED_ROOT=/usr/local/other/Fortuna-2_5.noHDF5 MPPDIR=/usr/local/other/Fortuna-2_5.noHDF5/Linux FFTW_ROOT=/discover/nobackup/mkelley5/fftw-3.2.2";

$extraFlags{intel} = "";
$extraFlags{gfortran} ="";
$extraFlags{nag} ="";

# -----------------------------------------------------------------------------
sub createTemporaryCopy 
{
  my $referenceDir = shift;
  my $tempDir = shift;
  my $branch = shift;
  my $commandString = "git clone -b $branch $referenceDir $tempDir";
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
  #print "2) COMPILER: $env->{COMPILER}\n";
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
  my $gitroot          = $env -> {GIT_REPOSITORY};
  my $branch           = $env -> {BRANCH};

  my $commandString = <<EOF;
  pushd $scratchDirectory
  git clone -b $branch $gitroot $branch
  popd
EOF
  print "gitCheckout: $commandString\n";
  return (CommandEntry -> new({COMMAND => $commandString}))
}
1;

