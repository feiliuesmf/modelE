use Env;
use CommandPool;
use CommandEntry;
use RegressionTools;

open(LOG,">nightlyTests.log", O_RDWR|O_CREAT, 0664);
# Make the LOG hot: <http://www.thinkingsecure.com/docs/tpj/issues/vol3_3/tpj0303-0002.html>
{
  my $ofh = select LOG;
  $| = 1;
  select $ofh;
}

$cfgFile = $ARGV[0];
print "Configuration file name: $cfgFile\n";

# Get rundecks, compiler, level settings from configuration file
my $localCfgFile = $ENV{MODELROOT} . "/exec/testing/" . "$cfgFile";
print "Local file name: $localCfgFile\n";
eval { require "$localCfgFile"};
if ($@) {
    print "Failed to load, because : $@"
}

# get references to arrays:
my $compilers = \@comps;

my $scratchDir = $ENV{REGSCRATCH}; 
my $resultsDir = $ENV{REGRESULTS};

my $rundecks;
my $branch;
($sec, $min, $hour, $day, $mon, $yrOffset, $dyOfWeek, $dayOfYr, $dyltSav) = localtime();

$branch = $gitbranch;
$rundecks = \@decks;
$gitdir = $ENV{MODELROOT};
print "Test $branch branch on directory $gitdir\n";
`cd $gitdir; git pull; cd -`;

my $reference = "$scratchDir/$branch";

my $env = {};
$env->{BRANCH} = $branch;
$env->{"intel"} = getEnvironment("intel", $scratchDir, $branch);
$env->{"gfortran"} = getEnvironment("gfortran", $scratchDir, $branch);

my $resolutions = {};
# HEAD rundecks
$resolutions->{EC12} = "8x10";
$resolutions->{EM20} = "4x5";
$resolutions->{E1oM20} = "4x5";
$resolutions->{E4F40} = "2x2.5";
$resolutions->{E4TcadF40} = "2x2.5";
$resolutions->{E4arobio_h4c} = "2x2.5";
$resolutions->{E4arobio_g6c} = "2x2.5";
$resolutions->{SCMSGPCONT} = "0"; # single column model 
$resolutions->{E4C90L40} = "CS"; # cubed sphere
# AR5 rundecks
$resolutions->{E_AR5_CADI} = "2x2.5";

# Save settings for diffreport
&saveForDiffreport($branch, $cfgFile);

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
#my $clean = CommandEntry->new({COMMAND => "rm -rf $scratchDir/* *.o[0-9]* $resultsDir/*/*;"});
my $clean = CommandEntry->new({COMMAND => "ls $scratchDir/* $resultsDir/*;"});
my $git = CommandEntry->new(gitCheckout($env->{"intel"})); # Compiler is not important here

$pool->add($clean);
$pool->add($git);

foreach $compiler (@$compilers) 
{
  $pool->add(writeModelErcFile($env->{$compiler}, $git));
  unless (-d $env->{$compiler}->{RESULTS_DIRECTORY}."/$compiler") 
  {
    mkdir "$env->{$compiler}->{RESULTS_DIRECTORY}/$compiler" or die $!;
  }
}

print "Loop over all configurations...\n";
foreach my $rundeck (@$rundecks) 
{ 
  foreach $compiler (@{$useCases->{$rundeck}->{COMPILERS}}) 
  {
    $env->{$compiler}->{RUNDECK} = $rundeck;
    foreach $configuration (@{$useCases->{$rundeck}->{CONFIGURATIONS}}) 
    {
      $env->{$compiler}->{CONFIGURATION} = $configuration;
      my $tempDir="$scratchDir/$compiler/$rundeck.$configuration";

      my $copy  = createTemporaryCopy($reference, $tempDir);
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
print LOG "Completed nightly regression tests.\n";
close(LOG);

