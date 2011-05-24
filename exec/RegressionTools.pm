#package RegressionTools;
use CommandEntry;
use Env;

my $HYCOM = "E1fzhyc";

my $extraFlags;
$extraFlags{SERIAL}   = "";
$extraFlags{MPI}      = "ESMF=YES NPES=\$npes";
$extraFlags{OPENMP}   = "EXTRA_FFLAGS=-mp MP=YES NPROCS=\$npes";
$extraFlags{SERIALMP} = "EXTRA_FFLAGS=-mp"; # Intel kludge for matching serial and OpenMP.

$extraFlags{E1M20} ="";
$extraFlags{E1F20} ="";
$extraFlags{E1oM20} ="";
$extraFlags{E001tr} ="";
$extraFlags{$HYCOM} ="EXTRA_FFLAGS+=-DCHECK_OCEAN_HYCOM";

$extraFlags{intel} = "";
$extraFlags{gfortran} ="";

sub createTemporaryCopy {
  my $referenceDir = shift;
  my $tempDir = shift;
  my $commandString = "/usr/local/other/git/1.7.3.4_GNU/bin/git clone $referenceDir $tempDir;";
  return (CommandEntry -> new({COMMAND => $commandString}));

  }

# Build rundeck and compile
# Do a vclean just in case receiving a dirty copy.
# Creat CMPE002P if SERIAL configuration - used for verification
sub compileRundeck {
  my $env = shift;
  my $installDir = shift;

  my $compiler = $env->{COMPILER};
  my $rundeck  = $env->{RUNDECK};
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

  my $commandString = <<EOF;
  export MODELERC=$MODELERC;
  echo "MODELERC is $MODELERC";
  echo "CONTENTS: ";
  cat $MODELERC;
  cd $installDir/decks;
  make rundeck RUN=$expName RUNSRC=$rundeck;
  make vclean RUN=$expName;
  make -j gcm RUN=$expName $flags COMPILER=$compiler;
EOF

  my $binDir = $expName . "_bin";
  print "Compile? $commandString \n";
  return (CommandEntry -> new({COMMAND => $commandString, QUEUE => "", STDOUT_LOG_FILE => "$logFile", COMPILER => $compiler, MODELERC=>$MODELERC }));
}

sub runConfiguration {
    my $env = shift;
    my $installDir = shift;
    my $npes = shift;

    my $compiler = $env->{COMPILER};
    my $rundeck    = $env->{RUNDECK};
    my $configuration = $env->{CONFIGURATION};
    my $resultsDir = $env->{RESULTS_DIRECTORY};
    $resultsDir .="/$compiler";

    print "results dir $resultsDir \n";

    my $flags = "$extraFlags{$configuration}";
    $flags =~ s/(\$npes)/$npes/eeg;
    my $expName = "$rundeck.$configuration.$compiler";

    my $suffix;
    my $mpiArgs;
    my $logFile = "$resultsDir/$expName.runlog";

    if ($configuration eq "SERIAL") {
	$suffix = "";
	$mpiArgs = "";
    }
    else {
	$suffix = ".np=$npes";
	$mpiArgs = "-np $npes";
    }
    
    my $MODELERC = $env->{MODELERC};

    my $run1hr = "make rundeck RUN=$expName RUNSRC=$rundeck OVERWRITE=YES; make setup RUN=$expName $flags; MODELERC=$MODELERC ../exec/runE $expName -np $npes -cold-restart";
    my $run1dy = "$installDir/exec/editRundeck.sh $expName 48 2 1; make setup RUN=$expName $flags; MODELERC=$MODELERC ../exec/runE $expName -np $npes -cold-restart";
    my $restart = "pushd $expName; cp fort.2.nc fort.1.nc ; ./$expName $mpiArgs  ; popd";

    if ($configuration eq "MPI" or $configuration eq "OPENMP") {$continue1Day = "./$expName -np $npes ";}

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
  return (CommandEntry -> new({COMMAND => $commandString, QUEUE => "", STDOUT_LOG_FILE => "$logFile", NUM_PROCS => $npes, COMPILER => $compiler} ));
}


sub writeModelErcFile {
    my $env = shift;
    my $modelerc = $env->{MODELERC};
    my $commandString .= "mkdir -p $env->{SCRATCH_DIRECTORY}/$env->{COMPILER}\n"; 
    $commandString .= "rm $modelerc\n";
    
    while (my ($var, $value) = each(%$env) ) {
        $commandString .= "echo $var=$value >> $modelerc\n";
    }
    $commandString .= "mkdir -p $env->{DECKS_REPOSITORY} $env->{CMRUNDIR} $env->{SAVEDISK} $env->{EXECDIR} \n";
    return (CommandEntry -> new({COMMAND => $commandString}))
}

sub gitCheckout {
  my $env = shift;

  my $scratchDirectory = $env -> {SCRATCH_DIRECTORY};
  my $gitroot          = $env -> {GITROOT};

  my $commandString = <<EOF;
(pushd $scratchDirectory
/usr/local/other/git/1.7.3.4_GNU/bin/git clone $gitroot
popd)
EOF
  return (CommandEntry -> new({COMMAND => $commandString}))
}

sub getEnvironment {
    my $compiler = shift;
    my $scratchDir = shift;

    if ($compiler eq "intel") {
	return getIntelEnvironment($scratchDir);
    }
    else {
	return getGfortranEnvironment($scratchDir);
    }
}

sub getIntelEnvironment{
    my $scratchDir = shift;

    my $env = {};

    $env->{SCRATCH_DIRECTORY}=$scratchDir;
    $env->{BASELINE_DIRECTORY}="$ENV{NOBACKUP}/modelE_baseline";
    $env->{RESULTS_DIRECTORY} = $ENV{NOBACKUP}."/regression_results";
#    $env->{GITROOT}="simplex.giss.nasa.gov:/giss/gitrepo/modelE.git";
    $env->{GITROOT}="/home/modele/modelE";
    $env->{DECKS_REPOSITORY}="$scratchDir/decks_repository";
    $env->{CMRUNDIR}="$scratchDir/cmrun";
    $env->{EXECDIR}="$scratchDir/exec";
    $env->{SAVEDISK}="$scratchDir/savedisk";
    $env->{GCMSEARCHPATH}="/discover/nobackup/projects/giss/prod_input_files";
    $env->{MP}="no";
    $env->{OVERWRITE}="YES";
    $env->{OUTPUT_TO_FILES}="YES";
    $env->{VERBOSE_OUTPUT}="YES";
    $env->{BASELIBDIR5}="/usr/local/other/esmf510/Linux";
    $env->{MPIDISTR}="intel";
    $env->{COMPILER}="intel";
    $env->{ESMF_BOPT}="O";
    $env->{NETCDFHOME}="/usr/local/other/netcdf/3.6.2_intel-11.0.083";
    $env->{PNETCDFHOME}="/usr/local/other/pnetcdf/intel11.1.072_impi3.2.2.006";
    $env->{MODELERC}="$scratchDir/intel/modelErc.intel";
    return $env;
}

sub getGfortranEnvironment {
    my $scratchDir = shift;

    my $env = {};
    
    $env->{SCRATCH_DIRECTORY}=$scratchDir;
    $env->{BASELINE_DIRECTORY}="$ENV{NOBACKUP}/modelE_baseline";
    $env->{RESULTS_DIRECTORY} = $ENV{NOBACKUP}."/regression_results";
#    $env->{GITROOT}="simplex.giss.nasa.gov:/giss/gitrepo/modelE.git";
    $env->{GITROOT}="/home/modele/modelE";
    $env->{DECKS_REPOSITORY}="$scratchDir/decks_repository";
    $env->{CMRUNDIR}="$scratchDir/cmrun";
    $env->{EXECDIR}="$scratchDir/exec";
    $env->{SAVEDISK}="$scratchDir/savedisk";
    $env->{GCMSEARCHPATH}="/discover/nobackup/projects/giss/prod_input_files";
    $env->{MP}="no";
    $env->{OVERWRITE}="YES";
    $env->{OUTPUT_TO_FILES}="YES";
    $env->{VERBOSE_OUTPUT}="YES";
    $env->{BASELIBDIR5}="/usr/local/other/esmf5/gcc4.5_openmpi-1.4.2/Linux";
    $env->{MPIDIR}="/usr/local/other/openMpi/gcc-4.5";
    $env->{MPIDISTR}="openmpi";
    $env->{COMPILER}="gfortran";
    $env->{ESMF_BOPT}="O";
    $env->{NETCDFHOME}="/usr/local/other/netcdf/3.6.2_gcc4.5";
    $env->{PNETCDFHOME}="/usr/local/other/pnetcdf/gcc4.5_openmpi-1.4.2";
    $env->{MODELERC}="$scratchDir/gfortran/modelErc.gfortran";
    return $env;
}

sub checkConsistency {
    my $env = shift;
    my $rundeck = shift;
    my $useCases = shift;

    my $results = {};
    $results->{ARE_CONSISTENT} = 1; # True unless proved otherwise
    $results->{ALL_COMPLETED} = 1; # True unless proved otherwise
    $results->{NEW_SERIAL} = 0; 
    $results->{MESSAGES} = "";
    $results->{RESTART_CHECK} = 1; # True unless proved otherwise

    $compiler = $env->{COMPILER};

    my $compare = "/discover/nobackup/projects/giss/exec/diffreport";

    foreach my $configuration (@{$useCases->{CONFIGURATIONS}}) {

	my $tempDir="$scratchDir/$compiler/$rundeck.$configuration.$compiler";
	    
	@peList = (1);
	if ($configuration eq "MPI" or $configuration eq "OPENMP") {
	    @peList = @{$useCases->{NUM_MPI_PROCESSES}};
	}

	my $resultsDir = $env->{RESULTS_DIRECTORY} . "/$compiler";
	    
	foreach my $npes (@peList) {
	    my $suffix;
	    if    ($configuration eq "MPI")    {$suffix = ".np=$npes";}
	    else {$suffix = "";}
	    
	    foreach my $duration (@{$useCases->{DURATIONS}}) {
		
		my $reference;
		if    ($configuration eq "MPI" or $duration eq "restart")    {$reference = "$rundeck.SERIAL";}
		else {$reference = "$env->{BASELINE_DIRECTORY}/$compiler/$rundeck.SERIAL";}
	    
		my $outfile = "$resultsDir/$rundeck.$configuration.$compiler.$duration$suffix";
		print LOG "Looking for $outfile \n";
		if (-e $outfile) {
		    my $referenceOutput;
		    if ($duration == "restart") {
			$referenceOutput = "$reference.$compiler.1dy";
		    }
		    else {
			$referenceOutput = "$reference.$compiler.$duration";
		    }
		    my $testOutput = "$rundeck.$configuration.$compiler.$duration$suffix";
		    my $numLinesFound = `cd $resultsDir; $compare $referenceOutput $testOutput is_npes_reproducible | wc -l`;
		    chomp($numLinesFound);
		    print LOG "CHOMP:: reference: $reference \n";
		    print LOG "CHOMP:: referenceOutput: $referenceOutput \n";
		    print LOG "CHOMP:: testOutput: $testOutput \n";
		    print LOG "CHOMP:: $numLinesFound \n";
		    if ($numLinesFound > 0) {
			if ($configuration eq "MPI") {
			    $results->{MESSAGES} .= "   FAILURE - Inconsistent results for rundeck $rundeck, compiler $compiler, and duration $duration on $npes npes.\n";
			    $results->{ARE_CONSISTENT} = 0; #failure
			}
			else {
			    if ($duration eq "restart") {
				$results->{MESSAGES} .= "   FAILURE - Serial restart is not consistent.\n";
				$results->{RESTART_CHECK} = 0;
			    }
			    else {
				$results->{MESSAGES} .= "   WARNING - Serial run of rundeck $rundeck using compiler $compiler is different than stored results for duration $duration. \n";
				$results->{NEW_SERIAL} = 1; # change detected
			    }
			}
		    }
		}
		else {
		    print LOG " but it was not found \n";
		    $results->{MESSAGES} .= "   FAILURE - expected output file does not exist: $outfile\n";
		    $results->{ALL_COMPLETED} = 0; # failure
		    $results->{ARE_CONSISTENT} = 0; # failure
		    $results->{RESTART_CHECK} = 0; # failure
		    last;
		}
	    }
	}
    }
    return $results;
}

1;
