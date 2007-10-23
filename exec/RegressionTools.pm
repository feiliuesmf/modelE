#package RegressionTools;
use CommandEntry;

my $standardFlags = "SETUP_FLAGS=-wait";

my $extraFlags;
$extraFlags{SERIAL}   = "";
$extraFlags{MPI}      = "ESMF=YES NPES=\$npes";
$extraFlags{OPENMP}   = "EXTRA_FFLAGS=-mp MP=YES NPROCS=\$npes";
$extraFlags{SERIALMP} = "EXTRA_FFLAGS=-mp"; # Intel kludge for matching serial and OpenMP.

my $numLinesCMP;
$numLinesCMP{E1M20}  = 120;
$numLinesCMP{E1F20}  = 120;
$numLinesCMP{E1oM20} = 110;
$numLinesCMP{E001tr} = 123;

sub setModuleEnvironment {
#    require perlModule;
    require "$ENV{MODULESHOME}/init/perl";
    module (purge);
#    module (load, "comp/intel-9.1.042", "mpi/scali-5.3", "lib/mkl-9.0.017");
    module (load, "comp/intel-9.1.049", "mpi/scali-5.3", "lib/mkl-9.0.017");
}

sub createTemporaryCopy {
  my $env = shift;
  my $tempDir = shift;

  my $referenceDir = $env -> {REFERENCE_DIRECTORY};
  my $commandString = "mkdir -p $tempDir;  cp -r -u $referenceDir/* $tempDir;";
  return (CommandEntry -> new({COMMAND => $commandString}));

  }

# Build rundeck and compile
# Do a vclean just in case receiving a dirty copy.
# Creat CMPE002P if SERIAL configuration - used for verification
sub compileRundeck {
  my $env = shift;
  my $installDir = shift;

  my $rundeck    = $env -> {RUNDECK};
  my $configuration = $env -> {CONFIGURATION};
  my $resultsDir = $env -> {RESULTS_DIRECTORY};

  my $flags = "$standardFlags $extraFlags{$configuration}";
  $flags =~ s/(\$npes)/1/eeg;

  my $expName;
  if (@_) {$expName = shift}
  else {$expName = "$rundeck.$configuration"};

  my $commandString = <<EOF;
  export MODELERC;
  cd $installDir/decks;
  make rundeck RUN=$expName RUNSRC=$rundeck;
  make vclean RUN=$expName;
  make gcm RUN=$expName $flags;
EOF

  my $binDir = $expName . "_bin";
  if ($configuration eq SERIAL) {
    $commandString .= "make aux RUN=$expName $flags;\n";
    $commandString .= "cp $binDir/CMPE002P $resultsDir/CMPE002P.$rundeck;\n";
    $commandString .= "cp $binDir/CMPE002 $resultsDir/CMPE002.$rundeck;\n";
  }

  return (CommandEntry -> new({COMMAND => $commandString, QUEUE => "datamove", STDOUT_LOG_FILE => "$resultsDir/$expName.buildlog"}));
}

sub runConfiguration {
    my $env = shift;
    my $installDir = shift;
    my $npes = shift;

    my $rundeck    = $env -> {RUNDECK};
    my $configuration = $env -> {CONFIGURATION};
    my $resultsDir = $env -> {RESULTS_DIRECTORY};

    my $flags = "$standardFlags $extraFlags{$configuration}";
    $flags =~ s/(\$npes)/$npes/eeg;
    my $expName = "$rundeck.$configuration";

    my $suffix;
    if ($configuration eq "SERIAL" or $configuration eq "SERIALMP") {
	$suffix = "";
    }
    else {
	$suffix = ".np=$npes";
    }
    
    my $continue1Day = "./$expName -r ";
    if ($configuration eq "MPI" or $configuration eq "OPENMP") {$continue1Day = "./$expName -np $npes -r ";}

    my $commandString = <<EOF;
    echo "Using flags: $flags "
    export MODELERC;
    cd $installDir/decks;
    make setup_nocomp RUN=$expName $flags;
    cd $expName;
    # Remove old result so that we notice if things fail really badly
    rm $resultsDir/$expName.1hr$suffix;
    cp fort.2 $resultsDir/$expName.1hr$suffix;
    $continue1Day;
    # Remove old result so that we notice if things fail really badly
    rm $resultsDir/$expName.1dy$suffix;
    cp fort.2 $resultsDir/$expName.1dy$suffix;
EOF
  return (CommandEntry -> new({COMMAND => $commandString, QUEUE => "", STDOUT_LOG_FILE => "$resultsDir/$expName.runlog", NUM_PROCS => $npes}));
}

sub runSetup {
  my $env = shift;

  my $installDir    = $env -> {INSTALL_DIR};
  my $rundeck       = $env -> {RUNDECK};
  my $configuration = $env -> {CONFIGURATION};
  my $npes          = $env -> {NUM_PROCS};
  my $resultsDir    = $env -> {RESULTS_DIR};

  my $expName;
  if (@_) {$expName = shift}
  else {$expName = "$rundeck.$configuration"};

  my $commandString = <<EOF;
cd $dir/decks;
make setup_nocompile RUN=$expName $standardFlags;
cp $expName/fort.2 $resultsDir/$expName.1hr;
pushd $expName;
./$expName -r;
cp $expName/fort.2 $resultsDir/$expName.1dy;
EOF

  return (CommandEntry -> new({COMMAND => $commandString}))
}

sub checkResults {
  my $env = shift;
  my $tempDir = shift;
  my $npes = shift;

  my $installDir    = $env -> {INSTALL_DIR};
  my $rundeck       = $env -> {RUNDECK};
  my $configuration = $env -> {CONFIGURATION};
  my $resultsDir    = $env -> {RESULTS_DIR};
  my $cmp           = "$resultsDir/CMPE002";
  my $expName;

  if (@_) {$expName = shift}
  else {$expName = "$rundeck.$configuration"};

  my $suffix;
  if ($configuration eq "SERIAL" or $configuration eq "SERIALMP") {
      $suffix = "";
  }
  else {
      $suffix = ".np=$npes";
  }

  my $commandString = <<EOF;
cd $resultsDir;
$cmp $rundeck.SERIAL.1hr $expName.1hr$suffix;
$cmp $rundeck.SERIAL.1dy $expName.1dy$suffix;
EOF

  return (CommandEntry -> new({COMMAND => $commandString, STDOUT_LOG_FILE => "$expName.log"}))
}

sub cvsCheckout {
#    my $what = shift;
  my $env = shift;

  my $scratchDirectory = $env -> {SCRATCH_DIRECTORY};
  my $cvsroot          = $env -> {CVSROOT};

  my $commandString = <<EOF;
(pushd $scratchDirectory
cvs -d $cvsroot co modelE
popd)
EOF

  return (CommandEntry -> new({COMMAND => $commandString}))
}

1;
