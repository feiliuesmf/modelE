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

@modelErcVariables = (DECKS_REPOSITORY, MP, CMRUNDIR, EXECDIR, OVERWRITE, 
		      OUTPUT_TO_FILES, VERBOSE_OUTPUT, SAVEDISK, GCMSEARCHPATH,
		      COMPILER, COMPILER_VERSION, BASELIBDIR5, ESMF_BOPT, MPIDISTR, NETCDFHOME);

my $env;
$env->{"intel"}->{SCRATCH_DIRECTORY}=$scratchDir;
$env->{"intel"}->{REFERENCE_DIRECTORY}="$scratchDir/modelE";
$env->{"intel"}->{BASELINE_DIRECTORY}="$ENV{NOBACKUP}/modelE_baseline";
$env->{"intel"}->{RESULTS_DIRECTORY} = $ENV{NOBACKUP}."/regression_results";
$env->{"intel"}->{GITROOT}="simplex.giss.nasa.gov:/giss/gitrepo/modelE.git";
$env->{"intel"}->{DECKS_REPOSITORY}="$scratchDir/decks_repository";
$env->{"intel"}->{CMRUNDIR}="$scratchDir/cmrun";
$env->{"intel"}->{EXECDIR}="$scratchDir/exec";
$env->{"intel"}->{SAVEDISK}="$scratchDir/savedisk";
$env->{"intel"}->{GCMSEARCHPATH}="/discover/nobackup/projects/giss/prod_input_files";
$env->{"intel"}->{MP}="no";
$env->{"intel"}->{OVERWRITE}="YES";
$env->{"intel"}->{OUTPUT_TO_FILES}="YES";
$env->{"intel"}->{VERBOSE_OUTPUT}="YES";
$env->{"intel"}->{BASELIBDIR5}="/usr/local/other/esmf510/Linux";
$env->{"intel"}->{MPIDISTR}="intel";
$env->{"intel"}->{COMPILER}="intel";
$env->{"intel"}->{ESMF_BOPT}="O";
$env->{"intel"}->{NETCDFHOME}="/usr/local/other/netcdf/3.6.2_intel-11.0.083";
$env->{"intel"}->{MODELERC}="$scratchDir/intel/modelErc.intel";

$env->{"gfortran"}->{SCRATCH_DIRECTORY}=$scratchDir;
$env->{"gfortran"}->{REFERENCE_DIRECTORY}="$scratchDir/modelE";
$env->{"gfortran"}->{BASELINE_DIRECTORY}="$ENV{NOBACKUP}/modelE_baseline";
$env->{"gfortran"}->{RESULTS_DIRECTORY} = $ENV{NOBACKUP}."/regression_results";
$env->{"gfortran"}->{GITROOT}="simplex.giss.nasa.gov:/giss/gitrepo/modelE.git";
$env->{"gfortran"}->{DECKS_REPOSITORY}="$scratchDir/decks_repository";
$env->{"gfortran"}->{CMRUNDIR}="$scratchDir/cmrun";
$env->{"gfortran"}->{EXECDIR}="$scratchDir/exec";
$env->{"gfortran"}->{SAVEDISK}="$scratchDir/savedisk";
$env->{"gfortran"}->{GCMSEARCHPATH}="/discover/nobackup/projects/giss/prod_input_files";
$env->{"gfortran"}->{MP}="no";
$env->{"gfortran"}->{OVERWRITE}="YES";
$env->{"gfortran"}->{OUTPUT_TO_FILES}="YES";
$env->{"gfortran"}->{VERBOSE_OUTPUT}="YES";
$env->{"gfortran"}->{BASELIBDIR5}="/usr/local/other/esmf5/gcc4.5_openmpi-1.4.2/Linux";
$env->{"gfortran"}->{MPIDIR}="/usr/local/other/openMpi/gcc-4.5";
$env->{"gfortran"}->{MPIDISTR}="openmpi";
$env->{"gfortran"}->{COMPILER}="gfortran";
$env->{"gfortran"}->{ESMF_BOPT}="O";
$env->{"gfortran"}->{NETCDFHOME}="/usr/local/other/netcdf/3.6.2_gcc4.5";
$env->{"gfortran"}->{MODELERC}="$scratchDir/gfortran/modelErc.gfortran";

my $rundecks = ["EM20"];
#my $rundecks = ["EM20", "E4F40", "E4TcadF40",
#		"E4arobio_h4c", "E4arobio_g6c", "SCMSGPCONT"];
my $compilers = ["intel", "gfortran"];
#my $compilers = ["intel"];

my $configurations;
$configurations->{"intel"}->{"EM20"}  = ["SERIAL", "MPI"];
$configurations->{"intel"}->{"E4F40"}  = ["SERIAL", "MPI"];
$configurations->{"intel"}->{"E4TcadF40"} = ["SERIAL", "MPI"];
$configurations->{"intel"}->{"E4arobio_h4c"} = ["SERIAL", "MPI"];
$configurations->{"intel"}->{"E4arobio_g6c"} = ["SERIAL", "MPI"];
$configurations->{"intel"}->{"SCMSGPCONT"} = ["SERIAL"];

$configurations->{"gfortran"}->{"EM20"}  = ["SERIAL", "MPI"];
$configurations->{"gfortran"}->{"E4F40"}  = ["SERIAL", "MPI"];
$configurations->{"gfortran"}->{"E4TcadF40"} = ["SERIAL", "MPI"];
$configurations->{"gfortran"}->{"E4arobio_h4c"} = ["SERIAL", "MPI"];
$configurations->{"gfortran"}->{"E4arobio_g6c"} = ["SERIAL", "MPI"];
$configurations->{"gfortran"}->{"SCMSGPCONT"}  = ["SERIAL"];

my $numProcessors;
    $numProcessors->{EM20} ->{OPENMP} = [1,4];
    $numProcessors->{E4F40} ->{OPENMP} = [1,4];
    $numProcessors->{SCMSGPCONT} ->{OPENMP} = [1,4];

$numProcessors->{SCMSGPCONT} ->{MPI}    = [1]; # no effective MPI for SCM case

my $level = "aggressive";

if ($level eq "gentle") { # 3 lats per proc
    $numProcessors->{EM20} ->{MPI}    = [1,4,15];
    $numProcessors->{E4F40} ->{MPI}    = [1,4,30];
    $numProcessors->{E4TcadF40} ->{MPI}    = [1,4,30];
    $numProcessors->{E4arobio_H4c} ->{MPI}    = [1,4,30];
    $numProcessors->{E4arobio_g6c} ->{MPI}    = [1,4,30];
}
elsif ($level eq "aggressive") { # aggressive - 2 lats
    $numProcessors->{EM20} ->{MPI}    = [1,4,23];
    $numProcessors->{E4F40} ->{MPI}    = [1,4,45];
    $numProcessors->{E4TcadF40} ->{MPI}    = [1,4,30];
    $numProcessors->{E4arobio_H4c} ->{MPI}    = [1,4,30];
    $numProcessors->{E4arobio_g6c} ->{MPI}    = [1,4,30];
}
else { # insane - 1 lat per proc
    $numProcessors->{EM20} ->{MPI}    = [1,44];
    $numProcessors->{E4F40} ->{MPI}    = [1,88];
    $numProcessors->{E4TcadF40} ->{MPI}    = [1,88];
    $numProcessors->{E4arobio_H4c} ->{MPI}    = [1,44];
    $numProcessors->{E4arobio_g6c} ->{MPI}    = [1,88];
}

my $pool = CommandPool->new();

my $git = CommandEntry->new(gitCheckout($env->{"intel"})); # Compiler is not important here
my $clean = CommandEntry->new({COMMAND => "rm -rf $scratchDir/* regression.o*;"});

$pool->add($clean);
$pool->add($git);

# This should only be necessary if compiler==gfortran
$ENV{PATH}="/usr/local/other/gcc/4.5/bin:/usr/local/other/openMpi/gcc-4.5/bin:".$ENV{PATH};
$ENV{LD_LIBRARY_PATH}="/usr/local/other/gcc/4.5/lib64:/usr/local/other/openMpi/gcc-4.5/lib:".$ENV{LD_LIBRARY_PATH};

foreach my $compiler (@$compilers) {
    print LOG "Using compiler $compiler\n";

    $pool->add(writeModelErcFile($env->{$compiler}, $clean));

    unless (-d $env->{$compiler}->{RESULTS_DIRECTORY}."/$compiler") {
	mkdir "$env->{$compiler}->{RESULTS_DIRECTORY}/$compiler" or die $!;
    }

    foreach my $rundeck (@$rundecks) {
	print LOG " ... processing rundeck: $rundeck \n";
	$env->{$compiler}->{RUNDECK}=$rundeck;
	foreach my $configuration (@{$configurations->{$compiler}->{$rundeck}}) {
	    $env->{$compiler}->{CONFIGURATION}=$configuration;
	    print LOG "       ... processing configuration $configuration\n";
	    my $tempDir="$scratchDir/$compiler/$rundeck.$configuration";

	    my $copy  = createTemporaryCopy($env->{$compiler}, $tempDir);
	    $copy->{STDOUT_LOG_FILE} = "$env->{$compiler}->{RESULTS_DIRECTORY}/$compiler/$rundeck.$configuration.$compiler.buildlog";
	    $pool->add($copy, $git);
	    my $build = compileRundeck($env->{$compiler}, $tempDir, $compiler);
	    $pool->add($build, $copy);
	    
	    my $previous = $build;
	    
	    my @peList = (1);
	    if ($configuration eq "MPI" or $configuration eq "OPENMP") {
		@peList = @{$numProcessors->{$rundeck}->{$configuration}};
	    }
	    foreach my $npes (@peList) {
		my $run = runConfiguration($env->{$compiler}, $tempDir, $npes, $compiler);
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


foreach my $compiler (@$compilers) {
    print LOG "Checking results for compiler $compiler: \n";
    my $resultsDir = $env->{$compiler}->{RESULTS_DIRECTORY};
    $resultsDir .= "/$compiler";
    foreach my $rundeck (@$rundecks) {
	next if $rundeck eq "SERIAL" or $rundeck eq "SERIALMP";
	
	my $cmp = "/discover/nobackup/projects/giss/exec/diffreport";
	print LOG "    ... rundeck $rundeck: \n";
	my $consistent = 1; # unless proved otherwise
	my $newSerial  = 0; 
	foreach my $configuration (@{$configurations->{$compiler}->{$rundeck}}) {
	    my $tempDir="$scratchDir/$compiler/$rundeck.$configuration.$compiler";
	    my $reference;
	    if    ($configuration eq "MPI")    {$reference = "$rundeck.SERIAL";}
	    elsif ($configuration eq "OPENMP") {$reference = "$rundeck.SERIALMP";}
	    elsif ($configuration eq "SERIAL") {$reference = "$env->{$compiler}->{BASELINE_DIRECTORY}/$compiler/$rundeck.SERIAL";}
	    else {next;}
	    
	    my $suffix;
	    
	    print LOG "         ... configuration: $configuration \n";
	    
	    @peList = (1);
	    if ($configuration eq "MPI" or $configuration eq "OPENMP") {
		@peList = @{$numProcessors->{$rundeck}->{$configuration}};
	    }
	    
	    foreach my $npes (@peList){
		if    ($configuration eq "MPI" or $configuration eq "OPENMP")    {$suffix = ".np=$npes";}
		else {$suffix = "";}
		
		print LOG "     ... NPES=$npes ";
		
		foreach my $duration ("1hr", "1dy") {
		    
		    print LOG "$duration: ";
		    
		    $numLinesFound = `cd $resultsDir; $cmp $reference.$compiler.$duration $rundeck.$configuration.$compiler.$duration$suffix is_npes_reproducible | wc -l`;
		    chomp($numLinesFound);
		    
		    if ($numLinesFound > 0){
			print REPORT "Rundeck $rundeck ($configuration, $duration, $compiler) differs on $npes processors.\n";
			print REPORT "   diffreport shows $numLinesFound lines of output.\n";
			print REPORT `cd $resultsDir; $cmp $reference.$compiler.$duration $rundeck.$configuration.$compiler.$duration$suffix is_npes_reproducible`;
			print REPORT "\n";
			if ($configuration eq "SERIAL") {$newSerial = 1;}
			else {$consistent = 0;}
		    }
		    else {
			$checklist->{$rundeck}->{"$configuration.$compiler.$duration$suffix"} = "succeeded";
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
		`cd $resultsDir; cp $rundeck.SERIAL.$compiler.1* $env->{$compiler}->{BASELINE_DIRECTORY}/$compiler/.; `;
	    }
	}
	else {
	    foreach my $configuration (@{$configurations->{$compiler}->{$rundeck}}) {
		#my $buildLogTail = `tail -20 $resultsDir/$rundeck.$configuration.$compiler.buildlog`;
		#print REPORT "\nTail of Build Log for $rundeck.$configuration:\n";
		#print REPORT "$buildLogTail\n";
	    }
	}
	
    }
}

close REPORT;
# Mail report to mailing list
`mail -s "discover results" giss-modelE-regression\@lists.nasa.gov < Report`;

print LOG "Completed nightly regression tests.\n";
close(LOG);

