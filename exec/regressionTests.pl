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

my $scratchDir = $ENV{NOBACKUP}."/regression_scratch";
my $reference = "$scratchDir/modelE";

my $env = {};
$env->{"intel"} = getEnvironment("intel",$scratchDir);
$env->{"gfortran"} = getEnvironment("gfortran",$scratchDir);


my $resolutions = {};
$resolutions->{EC12} = "8x10";
$resolutions->{EM20} = "4x5";
$resolutions->{EF40} = "2x2.5";
$resolutions->{E4TcadF40} = "2x2.5";
$resolutions->{E4arobio_h4c} = "2x2.5";
$resolutions->{E4arobio_g6c} = "2x2.5";
$resolutions->{SCMSGPCONT} = "0"; #hack - serial only


my $rundecks = ["EM20", "E4F40", "E4TcadF40",
		"E4arobio_h4c", "E4arobio_g6c", "SCMSGPCONT"];
my $compilers = ["intel", "gfortran"];

#my $rundecks = ["EM20"];
#my $compilers = ["intel"];

#set defaults

my $level = "AGGRESSIVE";
my $numProcesses = {};

$numProcesses->{"0"}->{GENTLE}     = [];
$numProcesses->{"0"}->{AGGRESSIVE} = [];
$numProcesses->{"0"}->{INSANE}     = [];

$numProcesses->{"8x10"}->{GENTLE}     = [1,4];
$numProcesses->{"8x10"}->{AGGRESSIVE} = [1,4,12];
$numProcesses->{"8x10"}->{INSANE}     = [1,4,22];

$numProcesses->{"4x5"}->{GENTLE}     = [1,4];
$numProcesses->{"4x5"}->{AGGRESSIVE} = [1,4,23];
$numProcesses->{"4x5"}->{INSANE}     = [1,4,23,44];

$numProcesses->{"2x2.5"}->{GENTLE}     = [1,8];
$numProcesses->{"2x2.5"}->{AGGRESSIVE} = [1,8,45];
$numProcesses->{"2x2.5"}->{INSANE}     = [1,8,45,88];

my $useCases = {};
foreach my $rundeck (@$rundecks) {
    $useCases->{$rundeck}->{COMPILERS} = $compilers;
    $useCases->{$rundeck}->{CONFIGURATIONS} = ["SERIAL","MPI"];
    $useCases->{$rundeck}->{NUM_MPI_PROCESSES} = $numProcesses->{$resolutions->{$rundeck}}->{$level};
    $useCases->{$rundeck}->{DURATIONS} = ["1hr","1dy"];
}

# Override anything else here
$useCases->{"SCMSGPCONT"}->{CONFIGURATIONS} = ["SERIAL"];


my $pool = CommandPool->new();

my $clean = CommandEntry->new({COMMAND => "rm -rf $scratchDir/* regression.o*;"});
my $git = CommandEntry->new(gitCheckout($env->{"intel"})); # Compiler is not important here

$pool->add($clean);
$pool->add($git);

# This should only be necessary if compiler==gfortran
$ENV{PATH}="/usr/local/other/gcc/4.5/bin:/usr/local/other/openMpi/gcc-4.5/bin:".$ENV{PATH};
$ENV{LD_LIBRARY_PATH}="/usr/local/other/gcc/4.5/lib64:/usr/local/other/openMpi/gcc-4.5/lib:".$ENV{LD_LIBRARY_PATH};

foreach $compiler (@$compilers) {
    $pool->add(writeModelErcFile($env->{$compiler}, $clean));
    unless (-d $env->{$compiler}->{RESULTS_DIRECTORY}."/$compiler") {
	mkdir "$env->{$compiler}->{RESULTS_DIRECTORY}/$compiler" or die $!;
    }
}

foreach my $rundeck (@$rundecks) { 
    print LOG "Using rundeck $rundeck $env->{$compiler}->{RUNDECK} \n";

    foreach $compiler (@{$useCases->{$rundeck}->{COMPILERS}}) {
	print LOG "   ... using compiler $compiler \n";
	$env->{$compiler}->{RUNDECK} = $rundeck;

	foreach $configuration (@{$useCases->{$rundeck}->{CONFIGURATIONS}}) {
	    print LOG "       ... using configuration $configuration \n";
	    $env->{$compiler}->{CONFIGURATION} = $configuration;

	    my $tempDir="$scratchDir/$compiler/$rundeck.$configuration";
	    my $copy  = createTemporaryCopy($reference, $tempDir);
	    $copy->{STDOUT_LOG_FILE} = "$env->{$compiler}->{RESULTS_DIRECTORY}/$compiler/$rundeck.$configuration.$compiler.buildlog";
	    $pool->add($copy, $git);

	    my $build = compileRundeck($env->{$compiler}, $tempDir);
	    $pool->add($build, $copy);
	    
	    my $previous = $build;
	    
	    my @peList = (1);
	    if ($configuration eq "MPI" or $configuration eq "OPENMP") {
		@peList = @{$useCases->{$rundeck}->{NUM_MPI_PROCESSES}};
	    }
	    foreach my $npes (@peList) {
		my $run = runConfiguration($env->{$compiler}, $tempDir, $npes);
		$pool->add($run, $previous, $compiler);
		$previous = $run;
	    }
	}
    }
}

print LOG "\n***************************\n";
print LOG "Starting regression tests.\n";
print LOG "***************************\n\n";

$pool->run(*LOG, true);

open(REPORT,">Report");
print "Completed: $pool->{notCompleted} \n";
if ($pool->{notCompleted}) {
    print REPORT "***************************************\n";
    print REPORT "Warning: time expired on driver script.\n";
    print REPORT "***************************************\n";
}

my $compare = "/discover/nobackup/projects/giss/exec/diffreport";

foreach my $rundeck (@$rundecks) { 
    foreach $compiler (@{$useCases->{$rundeck}->{COMPILERS}}) {
	my $resultsDir = $env->{$compiler}->{RESULTS_DIRECTORY} . "/$compiler";

	my $consistent = 1; # unless proved otherwise
	my $newSerial  = 0; 

	foreach $configuration (@{$useCases->{$rundeck}->{CONFIGURATIONS}}) {
	    print LOG "         ... configuration: $rundeck.$compiler.$configuration \n";

	    my $tempDir="$scratchDir/$compiler/$rundeck.$configuration.$compiler";
	    my $reference;
	    if    ($configuration eq "MPI")    {$reference = "$rundeck.SERIAL";}
	    elsif ($configuration eq "SERIAL") {$reference = "$env->{$compiler}->{BASELINE_DIRECTORY}/$compiler/$rundeck.SERIAL";}
	    else {next;}
	    
	    my $suffix;
	    
	    @peList = (1);
	    if ($configuration eq "MPI" or $configuration eq "OPENMP") {
		@peList = @{$useCases->{$rundeck}->{NUM_MPI_PROCESSES}};
	    }
	    
	    foreach my $npes (@peList) {
		if    ($configuration eq "MPI")    {$suffix = ".np=$npes";}
		else {$suffix = "";}
		
		print LOG "     ... NPES=$npes ";
		
		foreach my $duration (@{$useCases->{$rundeck}->{DURATIONS}}) {

		    my $outfile = "$resultsDir/$rundeck.$configuration.$compiler.$duration$suffix";
		    print LOG "Looking for $outfile \n";
		    if (-e $outfile) {
			print LOG "  found it\n";
			my $numLinesFound = `cd $resultsDir; $compare $reference.$compiler.$duration $rundeck.$configuration.$compiler.$duration$suffix is_npes_reproducible | wc -l`;
			chomp($numLinesFound);
			if ($numLinesFound > 0) {
			    if ($configuration eq "MPI") {
				print REPORT "Inconsistent results for rundeck $rundeck, compiler $compiler, and duration $duration on $npes npes.\n";
				$consistent = 0; #failure
			    }
			    else {
				print REPORT "Serial run of rundeck $rundeck using compiler $compiler is different than stored results for duration $duration. \n";
				$newSerial = 1;
			    }
			}
		    }
		    else {
			print LOG " but it was not found \n";
			print REPORT "$rundeck.$configuration.$compiler.$duration$suffix does not exist.\n";
			$consistent = 0;
			last;
		    }
		}

	    }

	}

	if ($consistent) {
	    print REPORT "Rundeck $rundeck with compiler $compiler is strongly reproducible.\n";
	    if ($newSerial) {
		print REPORT "  ... However serial results for rundeck $rundeck have changed.  Assuming change is intentional.\n";
	    }
	}
	print LOG "\n";
    }
}

close REPORT;
# Mail report to mailing list
`mail -s "discover results" giss-modelE-regression\@lists.nasa.gov < Report`;

print LOG "Completed nightly regression tests.\n";
close(LOG);

