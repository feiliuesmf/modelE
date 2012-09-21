use Env;
use CommandPool;
use CommandEntry;
use RegressionTools;
use RegressionEnvs;

# This is the main modelE regression tests script. Its role is to setup the testing
# environment for rundecks specified in the configuration file (required argument).

$num_args = $#ARGV + 1;
if ($num_args != 1) {
  print "\nUsage: regressionTests.pl config_filename\n";
  exit;
}

$cfgFile = $ARGV[0];
print "Configuration file name: $cfgFile\n";

# This is the LOG file.
open(LOG,">modelETests.log", O_RDWR|O_CREAT, 0664); 

# Get rundecks, compiler, level settings from configuration file
eval { require "$cfgFile"};
if ($@) {
    print "Failed to load, because : $@"
}
# get references to configuration file options:
my $compilers = \@comps;
my $rundecks = \@decks;
my $branch = $gitbranch;
my $cleanScratch = $doCleanScratch;

# get ENVironmental variables and initialize other settings for this run
my $env = {};

$env->{BRANCH} = $branch;
$env = setupENVvariables($env);

foreach $compiler (@$compilers) 
{
   $env->{$compiler} = getEnvironment($env, $compiler, $branch);
}

# Save settings for diffreport
&saveForDiffreport($env, $cfgFile);

my $resolutions = {};
# HEAD rundecks
$resolutions->{E_AR5_C12} = "8x10";
$resolutions->{E4TcadC12} = "8x10";
$resolutions->{EM20} = "4x5";
$resolutions->{E1oM20} = "4x5";
$resolutions->{E4F40} = "2x2.5";
$resolutions->{E4TcadF40} = "2x2.5";
$resolutions->{E4TcadiF40} = "2x2.5";
$resolutions->{E4arobio_h4c} = "2x2.5";
$resolutions->{E4arobio_g6c} = "2x2.5";
$resolutions->{SCMSGPCONT} = "0"; # single column model 
$resolutions->{E4C90L40} = "CS"; # cubed sphere
# AR5 rundecks
$resolutions->{E_AR5_CADI} = "2x2.5";

my $numProcesses = {};

$numProcesses->{"0"}->{GENTLE}     = [];
$numProcesses->{"0"}->{AGGRESSIVE} = [];
$numProcesses->{"0"}->{INSANE}     = [];

$numProcesses->{"CS"}->{GENTLE}     = [6];
$numProcesses->{"CS"}->{AGGRESSIVE} = [6,48];
$numProcesses->{"CS"}->{INSANE}     = [6,84];

$numProcesses->{"8x10"}->{GENTLE}     = [4];
$numProcesses->{"8x10"}->{AGGRESSIVE} = [1,4,12];
$numProcesses->{"8x10"}->{INSANE}     = [1,4,22];

$numProcesses->{"4x5"}->{GENTLE}     = [4];
$numProcesses->{"4x5"}->{AGGRESSIVE} = [1,4,23];
$numProcesses->{"4x5"}->{INSANE}     = [1,4,23,44];

$numProcesses->{"2x2.5"}->{GENTLE}     = [8];
$numProcesses->{"2x2.5"}->{AGGRESSIVE} = [1,45];
$numProcesses->{"2x2.5"}->{INSANE}     = [1,8,45,88];

my $useCases = {};
foreach my $rundeck (@$rundecks) 
{
  $useCases->{$rundeck}->{COMPILERS} = $compilers;
  if ($rundeck =~ m/C90/ || $rundeck =~ m/AR5/) 
  {
    $useCases->{$rundeck}->{CONFIGURATIONS} = ["MPI"];
  }
  else 
  {
    $useCases->{$rundeck}->{CONFIGURATIONS} = ["SERIAL","MPI"];
  }
  $useCases->{$rundeck}->{NUM_MPI_PROCESSES} = $numProcesses->{$resolutions->{$rundeck}}->{$level};
  $useCases->{$rundeck}->{DURATIONS} = ["1hr","1dy","restart"];
}

# Override anything else here
$useCases->{"SCMSGPCONT"}->{CONFIGURATIONS} = ["SERIAL"];

foreach my $rundeck (@$rundecks) 
{ 
  foreach $compiler (@{$useCases->{$rundeck}->{COMPILERS}}) 
  {
    $env->{$compiler}->{BRANCH} = $branch;
  }
}

# Create a new command pool:
my $pool = CommandPool->new();

# If we are running ALL the tests at the same time we probably want to clean up 
# the scratch space
if ($cleanScratch eq 'YES') {
  my $clean = CommandEntry->new({COMMAND => "rm -rf *.o[0-9]* *.diff $env->{SCRATCH_DIRECTORY}/* $env->{RESULTS_DIRECTORY}/*/*;"});
  $pool->add($clean);
}
my $git = CommandEntry->new(gitCheckout($env)); # Compiler is not important here
$pool->add($git);


foreach $compiler (@$compilers) 
{
  $pool->add(writeModelErcFile($env->{$compiler}));
  unless (-d $env->{$compiler}->{RESULTS_DIRECTORY}."/$compiler") 
  {
    mkdir "$env->{$compiler}->{RESULTS_DIRECTORY}/$compiler" or die $!; 
  }
}

print "Loop over all configurations...\n";
my $reference = "$env->{SCRATCH_DIRECTORY}" . "/$branch";

foreach my $rundeck (@$rundecks) 
{ 
  foreach $compiler (@{$useCases->{$rundeck}->{COMPILERS}}) 
  {
    if ($compiler eq 'nag') {
      $ENV{PATH}="/discover/nobackup/ccruz/Baselibs/mvapich2_1.8/nag-5.3-886/lib:".$ENV{PATH};
      $ENV{LD_LIBRARY_PATH}="/discover/nobackup/ccruz/Baselibs/mvapich2_1.8/nag-5.3-886/lib:".$ENV{LD_LIBRARY_PATH};
    }
    $env->{$compiler}->{RUNDECK} = $rundeck;
    foreach $configuration (@{$useCases->{$rundeck}->{CONFIGURATIONS}}) 
    {
      $env->{$compiler}->{CONFIGURATION} = $configuration;
      my $tempDir="$env->{SCRATCH_DIRECTORY}/$compiler/$rundeck.$configuration";

      my $copy  = createTemporaryCopy($reference, $tempDir, $branch);
      $copy->{STDOUT_LOG_FILE} = "$env->{$compiler}->{RESULTS_DIRECTORY}/$compiler/$rundeck.$configuration.$compiler.buildlog";
      $pool->add($copy, $git);

      my $build = compileRundeck($env->{$compiler}, $tempDir);

      $pool->add($build, $copy);
	    
      my $previous = $build; # for dependencies
	    
      my @peList = (1);
      if ($configuration eq "MPI") 
      {
	@peList = @{$useCases->{$rundeck}->{NUM_MPI_PROCESSES}};
      }
      foreach my $npes (@peList) 
      {
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
print LOG "Completed modelE regression tests.\n";
close(LOG);

