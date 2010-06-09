use Env;
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
$env{BASELINE_DIRECTORY}="$ENV{NOBACKUP}/modelE_baseline";
$env{RESULTS_DIRECTORY} = $ENV{NOBACKUP}."/regression_results";

@modelErcVariables = (DECKS_REPOSITORY, MP, CMRUNDIR, EXECDIR, OVERWRITE, 
		      OUTPUT_TO_FILES, VERBOSE_OUTPUT, SAVEDISK, GCMSEARCHPATH,
		      COMPILER, COMPILER_VERSION, BASELIBDIR, ESMF_BOPT, MPIDISTR, NETCDFHOME);

$env{DECKS_REPOSITORY}="$scratchDir/decks_repository";
$env{MP}="no";
$env{CMRUNDIR}="$scratchDir/cmrun";
$env{EXECDIR}="$scratchDir/exec";
$env{OVERWRITE}="YES";
$env{OUTPUT_TO_FILES}="YES";
$env{VERBOSE_OUTPUT}="YES";
$env{SAVEDISK}="$scratchDir/savedisk";
$env{GCMSEARCHPATH}="/discover/nobackup/projects/giss/prod_input_files";
$env{BASELIBDIR}="/usr/local/other/baselibs/ESMF222rp3_NetCDF362b6_10.1.017_intelmpi/Linux";
$env{MPIDISTR}="intel";
$env{COMPILER}="intel";
$env{COMPILER_VERSION}="10";
$env{ESMF_BOPT}="O";
$env{NETCDFHOME}="/usr/local/other/netcdf/3.6.2_intel-10.1.013";

$ENV{MODELERC}="$scratchDir/modelErc.regression";
open (MODELERC, ">$ENV{MODELERC}") or die('could not create $ENV{MODELERC} file.\n');
foreach my $var (@modelErcVariables) {
    print MODELERC "$var=$env{$var}\n";
}
close MODELERC;

`mkdir -p $env{DECKS_REPOSITORY} $env{CMRUNDIR} $env{SAVEDISK} $env{EXECDIR}`;

my $rundecks = ["E1M20","E1oM20","E1F20","E001tr","E4F40", 
                 "test_E4TcadF40", "E4arobio_h4c", "E4arobio_g6c"];
my $compilers = ["intel", "gfortran"];

my $configurations;
$configurations -> {"intel"} -> {"E1M20"}  = ["SERIAL", "MPI"];
$configurations -> {"intel"} -> {"E1oM20"} = ["SERIAL", "MPI"];
$configurations -> {"intel"} -> {"E1F20"}  = ["SERIAL", "MPI"];
$configurations -> {"intel"} -> {"E4F40"}  = ["SERIAL", "MPI"];
$configurations -> {"intel"} -> {"E001tr"} = ["SERIAL", "MPI"];
$configurations -> {"intel"} -> {"test_E4TcadF40"} = ["SERIAL", "MPI"];
$configurations -> {"intel"} -> {"E4arobio_h4c"} = ["SERIAL", "MPI"];
$configurations -> {"intel"} -> {"E4arobio_g6c"} = ["SERIAL", "MPI"];

$configurations -> {"gfortran"} -> {"E1M20"}  = ["SERIAL"];
$configurations -> {"gfortran"} -> {"E1oM20"} = ["SERIAL"];
$configurations -> {"gfortran"} -> {"E1F20"}  = ["SERIAL"];
$configurations -> {"gfortran"} -> {"E4F40"}  = ["SERIAL"];
$configurations -> {"gfortran"} -> {"E001tr"} = ["SERIAL"];
$configurations -> {"gfortran"} -> {"test_E4TcadF40"} = [];
$configurations -> {"gfortran"} -> {"E4arobio_h4c"} = [];
$configurations -> {"gfortran"} -> {"E4arobio_g6c"} = [];

my $numProcessors;
    $numProcessors -> {E1M20}  -> {OPENMP} = [1,4];
    $numProcessors -> {E1oM20} -> {OPENMP} = [1,4];
    $numProcessors -> {E1F20}  -> {OPENMP} = [1,4];
    $numProcessors -> {E001tr} -> {OPENMP} = [1,4];
    $numProcessors -> {E4F40}  -> {OPENMP} = [1,4];

my $level = "insane";

if ($level eq "gentle") { # 3 lats per proc
    $numProcessors -> {E1M20}  -> {MPI}    = [1,4,15];
    $numProcessors -> {E1oM20} -> {MPI}    = [1,4,15];
    $numProcessors -> {E001tr} -> {MPI}    = [1,4,15];
    $numProcessors -> {E1F20}  -> {MPI}    = [1,4,30];
    $numProcessors -> {E4F40}  -> {MPI}    = [1,4,30];
    $numProcessors -> {test_E4TcadF40}  -> {MPI}    = [1,4,30];
    $numProcessors -> {E4arobio_H4c}  -> {MPI}    = [1,4,30];
    $numProcessors -> {E4arobio_g6c}  -> {MPI}    = [1,4,30];
}
elsif ($level eq "aggressive") { # aggressive - 2 lats
    $numProcessors -> {E1M20}  -> {MPI}    = [1,4,23];
    $numProcessors -> {E1oM20} -> {MPI}    = [1,4,23];
    $numProcessors -> {E001tr} -> {MPI}    = [1,4,23];
    $numProcessors -> {E1F20}  -> {MPI}    = [1,4,45];
    $numProcessors -> {E4F40}  -> {MPI}    = [1,4,45];
    $numProcessors -> {test_E4TcadF40}  -> {MPI}    = [1,4,30];
    $numProcessors -> {E4arobio_H4c}  -> {MPI}    = [1,4,30];
    $numProcessors -> {E4arobio_g6c}  -> {MPI}    = [1,4,30];
}
else { # insane - 1 lat per proc
    $numProcessors -> {E1M20}  -> {MPI}    = [1,44];
    $numProcessors -> {E1oM20} -> {MPI}    = [1,44];
    $numProcessors -> {E001tr} -> {MPI}    = [1,44];
    $numProcessors -> {E1F20}  -> {MPI}    = [1,88];
    $numProcessors -> {E4F40}  -> {MPI}    = [1,88];
    $numProcessors -> {test_E4TcadF40}  -> {MPI}    = [1,88];
    $numProcessors -> {E4arobio_H4c}  -> {MPI}    = [1,44];
    $numProcessors -> {E4arobio_g6c}  -> {MPI}    = [1,88];
}

setModuleEnvironment();
my $cvs = CommandEntry -> new(cvsCheckout(\%env));
my $clean = CommandEntry -> new({COMMAND => "rm -rf $scratchDir/*;"});
#my $clean = CommandEntry -> new({COMMAND => "ls;"});
my $pool = CommandPool->new();

$pool -> add($cvs);
$pool -> setFinal($clean);

$ENV{PATH}="/usr/local/other/gcc/4.5/bin:".$ENV{PATH};
$ENV{LD_LIBRARY_PATH}="/usr/local/other/gcc/4.5/lib64:".$ENV{LD_LIBRARY_PATH};

foreach my $compiler (@$compilers) {
    print LOG "Using compiler $compiler\n";
    print "$scratchDir/$compiler \n";
    unless (-d "$scratchDir/$compiler") {
	mkdir "$scratchDir/$compiler" or die $!;
    }
    unless (-d "$env{RESULTS_DIRECTORY}/$compiler") {
	mkdir "$env{RESULTS_DIRECTORY}/$compiler" or die $!;
    }
    foreach my $rundeck (@$rundecks) {
	print LOG " ... processing rundeck: $rundeck \n";
	$env{RUNDECK}=$rundeck;
	
	foreach my $configuration (@{$configurations -> {$compiler} -> {$rundeck}}) {
	    $env{CONFIGURATION}=$configuration;
	    print LOG "       ... processing configuration $configuration\n";
	    my $tempDir="$scratchDir/$compiler/$rundeck.$configuration";

	    my $copy  = createTemporaryCopy(\%env, $tempDir);
	    $copy -> {STDOUT_LOG_FILE} = "$env{RESULTS_DIRECTORY}/$compiler/$rundeck.$configuration.buildlog";
	    $pool -> add($copy, $cvs);
	    my $build = compileRundeck(\%env, $tempDir, $compiler);
	    $pool -> add($build, $copy);
	    
	    my $previous = $build;
	    
	    my @peList = (1);
	    if ($configuration eq "MPI" or $configuration eq "OPENMP") {
		@peList = @{$numProcessors -> {$rundeck} -> {$configuration}};
	    }
	    foreach my $npes (@peList) {
		my $run = runConfiguration(\%env, $tempDir, $npes, $compiler);
		$pool -> add($run, $previous);
		$previous = $run;
	    }
	}
    }
}

print LOG "\n***************************\n";
print LOG "Starting regression tests.\n";
print LOG "***************************\n\n";

$pool -> run(LOG, true);

open(REPORT,">Report");


my $LINES_EXPECTED;
#$LINES_EXPECTED -> {E1M20}   = [94,96];
#$LINES_EXPECTED -> {E1oM20}  = [120,122];
#$LINES_EXPECTED -> {E1F20}   = [94,96];
#$LINES_EXPECTED -> {E001tr}  = [102,104];
#$LINES_EXPECTED -> {$HYCOM}  = [139,141];

$LINES_EXPECTED -> {E1M20}   = [4,5];
$LINES_EXPECTED -> {E4F40}   = [4,15];
$LINES_EXPECTED -> {test_E4TcadF40} = [4,15];
$LINES_EXPECTED -> {E1oM20}  = [4,5];
$LINES_EXPECTED -> {E1F20}   = [4,5];
$LINES_EXPECTED -> {E001tr}  = [4,5];
$LINES_EXPECTED -> {$HYCOM}  = [10,12];
$LINES_EXPECTED -> {E1arobio_mpi}  = [10,12];

foreach my $compiler (@$compilers) {
    print LOG "Checking results for compiler $compiler: \n";
    my $resultsDir = $env{RESULTS_DIRECTORY};
    $resultsDir .= "/$compiler";
    foreach my $rundeck (@$rundecks) {
	next if $rundeck eq "SERIAL" or $rundeck eq "SERIALMP";
	
	my $cmp = "CMPE002P.$rundeck";
	print LOG "    ... rundeck $rundeck: \n";
	my $consistent = 1; # unless proved otherwise
	my $newSerial  = 0; 
	foreach my $configuration (@{$configurations -> {$compiler} -> {$rundeck}}) {
	    my $tempDir="$scratchDir/$compiler/$rundeck.$configuration.$compiler";
	    my $reference;
	    if    ($configuration eq "MPI")    {$reference = "$rundeck.SERIAL";}
	    elsif ($configuration eq "OPENMP") {$reference = "$rundeck.SERIALMP";}
	    elsif ($configuration eq "SERIAL") {$reference = "$env{BASELINE_DIRECTORY}/$compiler/$rundeck.SERIAL";}
	    else {next;}
	    
	    my $suffix;
	    
	    print LOG "         ... configuration: $configuration \n";
	    
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
		    
		    $numLinesFound = `cd $resultsDir; $cmp $reference.$duration $rundeck.$configuration.$duration$suffix | wc -l`;
		    chomp($numLinesFound);
		    
		    my $lineRange = $LINES_EXPECTED -> {$rundeck};
		    my $lineMin = $lineRange -> [0];
		    my $lineMax = $lineRange -> [1];
		    if (($numLinesFound <  ($lineRange -> [0])) or ($numLinesFound > ($lineRange -> [1]))) {
			print REPORT "Rundeck $rundeck ($configuration, $duration, $compiler) failed on $npes processors.\n";
			print REPORT "   CMPE002 (found $lineMin < $numLinesFound < $lineMax) \n";
#		    my $mess = `cd $resultsDir; $cmp $reference.$duration $rundeck.$configuration.$duration$suffix | head -100`;
#		    print REPORT "$mess\n";
			
			if ($configuration eq "SERIAL") {$newSerial = 1;}
			else {$consistent = 0;}
		    }
		    else {
			$checklist -> {$rundeck} -> {"$configuration.$duration$suffix"} = "succeeded";
			print LOG "succeeded ";
		    }
		}
		print LOG "\n";
	    }
	    print LOG "\n";
	}
	if ($consistent) {
	    print REPORT "Rundeck $rundeck with compiler $compiler is strongly reproducible.\n";
	    if ($newSerial) {
		print REPORT "  ... However serial results for rundeck $rundeck have changed.  Assuming change is intentional.\n";
		`cd $resultsDir; cp $rundeck.SERIAL.1* $env{BASELINE_DIRECTORY}/$compiler/.; `;
	    }
	}
	else {
	    foreach my $configuration (@{$configurations -> {$compiler} -> {$rundeck}}) {
		my $buildLogTail = `tail -20 $resultsDir/$rundeck.$configuration.buildlog`;
		print REPORT "\nTail of Build Log for $rundeck.$configuration:\n";
		print REPORT "$buildLogTail\n";
	    }
	}
	
    }
}

close REPORT;
# Mail report to mailing list
`mail -s "discover results" giss-modelE-regression\@lists.nasa.gov < Report`;

print LOG "Completed nightly regression tests.\n";
close(LOG);

