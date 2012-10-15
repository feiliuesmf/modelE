use CommandEntry;
use Env;

my $extraFlags;
$extraFlags{SERIAL}   = "";
$extraFlags{MPI}      = "MPI=YES";
$extraFlags{OPENMP}   = "EXTRA_FFLAGS=-mp MP=YES NPROCS=\$npes";
$extraFlags{SERIALMP} = "EXTRA_FFLAGS=-mp"; # Intel kludge for matching serial and OpenMP.

$extraFlags{E4TcadC12} ="";
$extraFlags{E_AR5_C12} ="";
$extraFlags{E4TcadF40} ="";
$extraFlags{E4TcadiF40} ="";
$extraFlags{E4arobio_g6c} ="";
$extraFlags{E4arobio_h4c} ="";
$extraFlags{EM20} ="";
$extraFlags{E4F40} ="";
$extraFlags{E1oM20} ="";
$extraFlags{E_AR5_CADI} ="ESMF=YES";
# Need to set CS variables in the configuration file:
$extraFlags{E4C90L40} ="ESMF=YES FVCUBED=YES FVCUBED_ROOT=/discover/nobackup/ccruz/GEOS5/Ganymed-1_0_BETA4.noHDF5 MPPDIR=/discover/nobackup/ccruz/GEOS5/Ganymed-1_0_BETA4.noHDF5/MPP/intel12_impi32.noHDF5 FFTW_ROOT=/discover/nobackup/mkelley5/fftw-3.2.2";

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
  my $configuration = $env -> {CONFIGURATION};

  my $expName;
  if (@_) {$expName = shift}
  else {$expName = "$rundeck.$configuration.$compiler"};

  my $branch   = $env->{BRANCH};
  my $debugFlags = $env -> {DEBUGFLAGS};

  my $resultsDir = $env -> {RESULTS_DIRECTORY};
  $resultsDir .="/$compiler";
  my $MODELERC = $env->{MODELERC};
 
  my $flags;
  if ($debugFlags eq 'N')
  {
    $flags = "$extraFlags{$configuration} $extraFlags{$rundeck} $extraFlags{$compiler}";
  }
  else
  {
    my $dFlags;
    if ($compiler eq 'intel') { $dFlags="\"-O0 -g -traceback\""; }
    elsif ($compiler eq 'gfortran') { $dFlags="\"-O0 -g -fbacktrace\""; }
    elsif ($compiler eq 'nag') { $dFlags="\"-O0 -g -gline\""; }
    $flags = "$extraFlags{$configuration} $extraFlags{$rundeck} $extraFlags{$compiler} EXTRA_FFLAGS=$dFlags";
  }
  $flags =~ s/(\$npes)/1/eeg;

  my $logFile = "$resultsDir/$expName.buildlog";
  unlink($logFile); # delete it

  my $commandString;
  if ($rundeck =~ m/^(E_AR5)/i)
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
    make clean RUN=$expName;
    make -j gcm RUN=$expName $flags COMPILER=$compiler;
EOF
  }


  my $binDir = $expName . "_bin";
  print "compileRundeck: $commandString \n";
  return (CommandEntry -> new({COMMAND => $commandString, QUEUE => "", 
          STDOUT_LOG_FILE => "$logFile", COMPILER => $compiler, 
          MODELERC=>$MODELERC, RUNDECK => $rundeck, BRANCH => $branch }));
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
  my $debugFlags = $env -> {DEBUGFLAGS};
  my $duration = $env -> {DURATION};
  my $resultsDir = $env->{RESULTS_DIRECTORY};
  $resultsDir .="/$compiler";
  my $expName = "$rundeck.$configuration.$compiler";
  my $MODELERC = $env->{MODELERC};

  print "BRANCH: $branch \n";
  print "DURATION: $duration \n";
  print "resultsDir: $resultsDir \n";

  my $flags;
  if ($debugFlags eq 'N')
  {
    $flags = "$extraFlags{$configuration} $extraFlags{$rundeck} $extraFlags{$compiler}";
  }
  else
  {
    my $dFlags;
    if ($compiler eq 'intel') { $dFlags="\"-O0 -g -traceback\""; }
    elsif ($compiler eq 'gfortran') { $dFlags="\"-O0 -g -fbacktrace\""; }
    elsif ($compiler eq 'nag') { $dFlags="\"-O0 -g -gline\""; }
    $flags = "$extraFlags{$configuration} $extraFlags{$rundeck} $extraFlags{$compiler} EXTRA_FFLAGS=$dFlags";
  }
  $flags =~ s/(\$npes)/1/eeg;

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
  my $run1hr;
  my $run1dy;
  my $restart;
  my $runCustom;
  if ($duration == 0) 
  {  
    $run1hr = createMakeCommand($env, $expName, $npes, $flags, 2);
    $run1dy = createMakeCommand($env, $expName, $npes, $flags, 48);
    $restart = "pushd $expName; cp fort.2.nc fort.1.nc ; ./$expName $mpiArgs  ; popd";
  }
  else
  {
    $runCustom = createMakeCommand($env, $expName, $npes, $flags, $duration,);
  }

  # Create a run command
  my $commandString;
  if ($duration == 0) 
  {  
    $commandString = <<EOF;
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
  }
  else
  {
    $commandString = <<EOF;
    export MODELERC=$MODELERC;
    cd $installDir/decks;
    rm -f $expName/fort.2.nc;
    $runCustom;
    if [ -e $expName/fort.2.nc ]; then
      cp $expName/fort.2.nc $resultsDir/$expName.$duration;
    else
      exit 1;
    fi
EOF
  }

  print "runConfiguration: $commandString\n";
  return (CommandEntry -> new({COMMAND => $commandString, QUEUE => "", 
          STDOUT_LOG_FILE => "$logFile", NUM_PROCS => $npes, 
          COMPILER => $compiler, RUNDECK => $rundeck , BRANCH => $branch }));
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

# -----------------------------------------------------------------------------
sub createMakeCommand
{
  my $env = shift;
  my $exp = shift;
  my $npes = shift;
  my $flags = shift;
  my $time = shift;

  my $deck = $env->{RUNDECK};
  my $rc = $env->{MODELERC};
  if ($time == 2 || $time == 48 ) 
  {
    if ($deck =~ m/^(E_AR5)/i)
    {
      return ("export MODELERC=$rc; 
             ../exec/editRundeck.sh $exp $time 2 1;
             make -f new.mk rundeck RUN=$exp RUNSRC=$deck OVERWRITE=YES; 
             make -f new.mk -j setup RUN=$exp $flags; 
             ../exec/runE_new $exp -np $npes -cold-restart; 
             rm -f $exp/run_status")
    }
    else
    {
      if ($time == 2) {
        return ("export MODELERC=$rc; 
           make rundeck RUN=$exp RUNSRC=$deck OVERWRITE=YES;  
           make -j setup RUN=$exp $flags; 
           ../exec/runE $exp -np $npes -cold-restart; 
           rm -f $exp/run_status")
      }
      else {
        return ("export MODELERC=$rc; 
           make rundeck RUN=$exp RUNSRC=$deck OVERWRITE=YES;  
           ../exec/editRundeck.sh $exp $time 2 1; 
           make -j setup RUN=$exp $flags; 
           ../exec/runE $exp -np $npes -cold-restart; 
           rm -f $exp/run_status")
      }
    }
  }
  else
  {
      my $date = $time / 48;
      return ("export MODELERC=$rc; 
           make rundeck RUN=$exp RUNSRC=$deck OVERWRITE=YES;  
           ../exec/editRundeck.sh $exp $time $date 0; 
           make -j setup RUN=$exp $flags; 
           ../exec/runE $exp -np $npes -cold-restart; 
           rm -f $exp/run_status")
  }

}

1;

