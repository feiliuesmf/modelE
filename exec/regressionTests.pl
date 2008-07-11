use CommandPool;
use CommandEntry;
use RegressionTools;

open(LOG,">nightlyTests.log");
# Make the LOG hot: <http://www.thinkingsecure.com/docs/tpj/issues/vol3_3/tpj0303-0002.html>
{
    my $ofh = select LOG;
    $| = 1;
    select $ofh;
}

$ENV{CVS_RSH}="ssh";
$env{CVSROOT}="simplex.giss.nasa.gov:/giss/cvsroot";
my $scratchDir = $ENV{NOBACKUP}."/regression_scratch";
$env{SCRATCH_DIRECTORY}=$scratchDir;
$env{REFERENCE_DIRECTORY}="$env{SCRATCH_DIRECTORY}/modelE";
$env{BASELINE_DIRECTORY}="$ENV{NOBACKUP}/regression_baseline";
$env{RESULTS_DIRECTORY} = $ENV{NOBACKUP}."/regression_results";

@modelErcVariables = (DECKS_REPOSITORY, MP, CMRUNDIR, EXECDIR, OVERWRITE, 
		      OUTPUT_TO_FILES, VERBOSE_OUTPUT, SAVEDISK, GCMSEARCHPATH,
		      COMPILER, ESMF_DIR, ESMF_BOPT, MPIDISTR, NETCDFHOME);

$env{DECKS_REPOSITORY}="$scratchDir/decks_repository";
$env{MP}="no";
$env{CMRUNDIR}="$scratchDir/cmrun";
$env{EXECDIR}="$scratchDir/exec";
$env{OVERWRITE}="YES";
$env{OUTPUT_TO_FILES}="YES";
$env{VERBOSE_OUTPUT}="YES";
$env{SAVEDISK}="$scratchDir/savedisk";
$env{GCMSEARCHPATH}="/discover/nobackup/projects/giss/prod_input_files";
$env{COMPILER}="Intel8";
$env{ESMF_DIR}="/discover/nobackup/projects/giss/esmf_2_2_ifort_9.1.042";
$env{MPIDISTR}="SCALI";
$env{ESMF_BOPT}="O";
$env{NETCDFHOME}="";

$ENV{MODELERC}="$scratchDir/modelErc.regression";
open (MODELERC, ">$ENV{MODELERC}") or die('could not create $ENV{MODELERC} file.\n');
foreach my $var (@modelErcVariables) {
    print MODELERC "$var=$env{$var}\n";
}
close MODELERC;

`mkdir -p $env{DECKS_REPOSITORY} $env{CMRUNDIR} $env{SAVEDISK} $env{EXECDIR}`;

my $HYCOM="E1fzhyc";
my $rundecks = ["E1M20","E1oM20","E1F20","E001tr",$HYCOM,"E1arobio_mpi"];
#my $rundecks = ["E1F20"];

my $configurations;

$configurations -> {"E1M20"}  = ["SERIAL", "SERIALMP", "OPENMP", "MPI"];
$configurations -> {"E1oM20"} = ["SERIAL", "SERIALMP", "OPENMP", "MPI"];
$configurations -> {"E1F20"}  = ["SERIAL", "SERIALMP", "OPENMP", "MPI"];
#$configurations -> {"E1F20"}  = ["SERIAL", "MPI"];
$configurations -> {"E001tr"} = ["SERIAL", "SERIALMP", "OPENMP", "MPI"];
$configurations -> {$HYCOM} = ["SERIAL", "SERIALMP", "OPENMP", "MPI"];

my $numProcessors;
    $numProcessors -> {E1M20}  -> {OPENMP} = [1,4];
    $numProcessors -> {E1oM20} -> {OPENMP} = [1,4];
    $numProcessors -> {E1F20}  -> {OPENMP} = [1,4];
    $numProcessors -> {E001tr} -> {OPENMP} = [1,4];
    $numProcessors -> {$HYCOM} -> {OPENMP} = [1,4];

my $level = "aggressive";

if ($level eq "gentle") { # 3 lats per proc
    $numProcessors -> {E1M20}  -> {MPI}    = [1,4,15];
    $numProcessors -> {E1oM20} -> {MPI}    = [1,4,15];
    $numProcessors -> {E001tr} -> {MPI}    = [1,4,15];
    $numProcessors -> {E1F20}  -> {MPI}    = [1,4,30];
    $numProcessors -> {$HYCOM}  -> {MPI}    = [1,2];
}
elsif ($level eq "aggressive") { # aggressive - 2 lats per proc
    $numProcessors -> {E1M20}  -> {MPI}    = [1,4,23];
    $numProcessors -> {E1oM20} -> {MPI}    = [1,4,15];
    $numProcessors -> {E001tr} -> {MPI}    = [1,4,23];
    $numProcessors -> {E1F20}  -> {MPI}    = [1,4,45];
    $numProcessors -> {$HYCOM}  -> {MPI}    = [1,2];
}
else { # insane - 1 lat per proc
    $numProcessors -> {E1M20}  -> {MPI}    = [1,4,46];
    $numProcessors -> {E1oM20} -> {MPI}    = [1,4,46];
    $numProcessors -> {E001tr} -> {MPI}    = [1,4,46];
    $numProcessors -> {E1F20}  -> {MPI}    = [1,4,90];
    $numProcessors -> {$HYCOM}  -> {MPI}    = [1];
}

setModuleEnvironment();
my $cvs = CommandEntry -> new(cvsCheckout(\%env));
my $clean = CommandEntry -> new({COMMAND => "rm -rf $scratchDir/*;"});
my $pool = CommandPool->new();

$pool -> add($cvs);
$pool -> setFinal($clean);

foreach my $rundeck (@$rundecks) {
    print LOG "Processing rundeck: $rundeck \n";
    $env{RUNDECK}=$rundeck;

    foreach my $configuration (@{$configurations->{$rundeck}}) {
	$env{CONFIGURATION}=$configuration;
	print LOG "... processing configuration $configuration\n";
	my $tempDir="$scratchDir/$rundeck.$configuration";
 
	my $copy  = createTemporaryCopy(\%env, $tempDir);
	$copy -> {STDOUT_LOG_FILE} = "$env{RESULTS_DIRECTORY}/rundeck.$configuration.buildlog";
	$pool -> add($copy, $cvs);

	my $build = compileRundeck(\%env, $tempDir);
	$pool -> add($build, $copy);

	my $previous = $build;

	my @peList = (1);
	if ($configuration eq "MPI" or $configuration eq "OPENMP") {
	    @peList = @{$numProcessors -> {$rundeck} -> {$configuration}};
	}

	foreach my $npes (@peList) {
	    my $run = runConfiguration(\%env, $tempDir, $npes);
	    $pool -> add($run, $previous);
	    $previous = $run;
	}
    }
}

print LOG "\n***************************\n";
print LOG "Starting regression tests.\n\n";

$pool -> run(LOG);

open(REPORT,">Report");


my $LINES_EXPECTED;
#$LINES_EXPECTED -> {E1M20}   = [94,96];
#$LINES_EXPECTED -> {E1oM20}  = [120,122];
#$LINES_EXPECTED -> {E1F20}   = [94,96];
#$LINES_EXPECTED -> {E001tr}  = [102,104];
#$LINES_EXPECTED -> {$HYCOM}  = [139,141];

$LINES_EXPECTED -> {E1M20}   = [4,5];
$LINES_EXPECTED -> {E1oM20}  = [4,5];
$LINES_EXPECTED -> {E1F20}   = [4,5];
$LINES_EXPECTED -> {E001tr}  = [4,5];
$LINES_EXPECTED -> {$HYCOM}  = [4,5];

foreach my $rundeck (@$rundecks) {
    next if $rundeck eq "SERIAL" or $rundeck eq "SERIALMP";

    my $cmp = "$env{RESULTS_DIRECTORY}/CMPE002P.$rundeck";
    print LOG "Checking results for rundeck $rundeck: \n";
    my $consistent = 1; # unless proved otherwise
    my $newSerial  = 0; 
    foreach my $configuration (@{$configurations -> {$rundeck}}) {
	my $tempDir="$scratchDir/$rundeck.$configuration";
	my $reference;
	if    ($configuration eq "MPI")    {$reference = "$rundeck.SERIAL";}
	elsif ($configuration eq "OPENMP") {$reference = "$rundeck.SERIALMP";}
	elsif ($configuration eq "SERIAL") {$reference = "$env{BASELINE_DIRECTORY}/$rundeck.SERIAL";}
	else {next;}

	my $suffix;

	print LOG "  ... configuration: $configuration \n";

	@peList = (1);
	if ($configuration eq "MPI" or $configuration eq "OPENMP") {
	    @peList = @{$numProcessors -> {$rundeck} -> {$configuration}};
	}

 	foreach my $npes (@peList){
	    if    ($configuration eq "MPI" or $configuration eq "OPENMP")    {$suffix = ".np=$npes";}
	    else {$suffix = "";}

	    print LOG "     ... NPES=$npes ";

	    foreach my $duration ("1hr", "1dy") {

		print LOG "$duration: ";
		
		$numLinesFound = `cd $env{RESULTS_DIRECTORY}; $cmp $reference.$duration $rundeck.$configuration.$duration$suffix | wc -l`;
		chomp($numLinesFound);

		my $lineRange = $LINES_EXPECTED -> {$rundeck};
		my $lineMin = $lineRange -> [0];
		my $lineMax = $lineRange -> [1];
		if (($numLinesFound <  ($lineRange -> [0])) or ($numLinesFound > ($lineRange -> [1]))) {
		    print REPORT "Rundeck $rundeck ($configuration, $duration) failed on $npes processors.\n";
		    print REPORT "   CMPE002 (found $lineMin < $numLinesFound < $lineMax) \n";
		    my $resultsDir = $env{RESULTS_DIRECTORY};
		    my $mess = `cd $resultsDir; $cmp $reference.$duration $rundeck.$configuration.$duration$suffix | head -100`;
		    print REPORT "$mess\n";

		    if ($configuration eq "SERIAL") {$newSerial = 1;}
		    else {$consistent = 0;}
		}
		else {
		    $checklist -> {$rundeck} -> {"$configuration.$durationr$suffix"} = "succeeded";
		    print LOG "succeeded ";
		}
	    }
	    print LOG "\n";
	}
	print LOG "\n";
    }
    if ($consistent) {
	print REPORT "Rundeck $rundeck is strongly reproducible.\n";
	if ($newSerial) {
	    print REPORT "  ... However serial results for rundeck $rundeck have changed.  Assuming change is intentional.\n";
	    `cd $env{RESULTS_DIRECTORY}; cp $rundeck.SERIAL.1* $env{BASELINE_DIRECTORY}/.; `;
	}
    }
    else {
	foreach my $configuration (@{$configurations -> {$rundeck}}) {
	    my $buildLogTail = `tail -20 $resultsDir/$rundeck.$configuration.buildlog`;
	    print REPORT "\nTail of Build Log for $rundeck.$configuration:\n";
	    print REPORT "$buildLogTail\n";
	}
    }

}

close REPORT;
# Mail report to mailing list
`mail -s "discover results" giss-modelE-regression\@sourcemotel.gsfc.nasa.gov < Report`;

print LOG "Completed nightly regresison tests.\n";
close(LOG);

