#!/usr/bin/perl
package CommandEntry;
use Env;
use Switch;

sub new 
{
  my $proto = shift;
  my ($args) = @_;

  my $class = ref($proto) || $proto;
  my $self  = {};
  
  $self -> {COMMAND} = " ";
  $self -> {NUM_PROCS} = 1;
  $self -> {QUEUE} = LOCAL; # or BATCH
  $self -> {DEPENDENCY} = ""; # place holder for list of references to other commands
  $self -> {STATUS} = PENDING;
  $self -> {STDOUT_LOG_FILE} = "LOG";
  $self -> {STDERR_LOG_FILE} = "ERR";
  $self -> {COMPILER} = "";
  $self -> {RUNDECK} = "";
  $self -> {BRANCH} = "";

  # override defaults with arguments
  while ( my ($key, $value) = each(%$args) ) 
  {
    $self -> {$key} = $value;
  }

  bless ($self, $class);
  return $self;
}

sub dependency 
{
  my $self       = shift;
  return ($self -> {DEPENDENCY});
}

sub status 
{
  my $self = shift;
  return ($self -> {STATUS});
}

sub isComplete 
{
  my $self       = shift;
  return ( $self -> {STATUS} eq COMPLETE );
}

sub isRunning 
{
  my $self       = shift;
  return ( $self -> {STATUS} eq RUNNING );
}

sub isPending 
{
  my $self       = shift;
  return ( $self -> {STATUS} eq PENDING );
}

sub isReady 
{
  my $self       = shift;  
  my $pending    = ($self -> isPending());
  my $dependency = ($self -> dependency());
  if ($dependency eq "") { return $pending }
  else {return ($pending and ($self -> dependency() -> isComplete()) )};
}

sub setStatus 
{
  my $self      = shift;
  my $newStatus = shift;
  $self -> {STATUS} = $newStatus;
}

sub setDependency 
{
  my $self = shift;
  my $newDependency = shift;
  $self -> {DEPENDENCY} = $newDependency;
}

sub runInBatch 
{
  my $self          = shift;
  my $commandString = shift;
  my $queue         = shift;
  my $queueString   = " ";

  if ($queue) 
  {
    $queueString = "-q $queue\n";
  }

  my $NODE_SIZE = 12;
  my $nCPUS = $self -> {NUM_PROCS};
  my $nodes = 1 + int(($nCPUS-1)/$NODE_SIZE);
  my $ncpus = $nCPUS;
  $ncpus = $NODE_SIZE if ($ncpus > $NODE_SIZE);


  my $compiler = $self->{COMPILER};
  my $branch   = $self->{BRANCH};
  my $jobname  = $self->{RUNDECK};

  switch ($jobname) 
  {
      case [ "EM20", "E1oM20", "E4F40", "SCMSGPCONT" ] 
      { $walltime = "1:00:00\n" }
      case [ "E4C90L40", "E_AR5_CADI" ] 
      { $walltime = "2:00:00\n" }
      case [ "E4TcadF40", "E4arobio_g6c", "E4arobio_h4c" ] 
      { $walltime = "3:00:00\n" }
      else
      { $walltime = "8:00:00\n"; }
  }

  print " runInBatch: COMPILER=$compiler, jobname=$jobname, BRANCH=$branch\n";

  my $script = <<EOF;
#!/bin/bash
#PBS -l select=$nodes:ncpus=12:mpiprocs=12
#PBS -l walltime=$walltime
#PBS -W group_list=s1001
#PBS -N $jobname
#PBS -j oe
#PBS $queueString
#PBS -V

cd \$PBS_O_WORKDIR

. /usr/share/modules/init/bash
module purge

EOF

  if ($branch =~ m/AR5/) 
  {
    if ($compiler eq "intel") 
    {

      $script .= <<EOF;
module load comp/intel-10.1.017 mpi/impi-3.2.2.006
EOF

    } 
    else 
    {
      $script .= <<EOF;
module load other/comp/gcc-4.5 other/mpi/openmpi/1.4.2-gcc-4.5
EOF
    }

  } 
  else 
  {
    if ($compiler eq "intel") 
    {
      $script .= <<EOF;
module load comp/intel-11.1.072 mpi/impi-3.2.2.006
EOF
    } 
    else 
    {
      $script .= <<EOF;
module load other/comp/gcc-4.6 other/mpi/mvapich2-1.4.1/gcc-4.6
EOF
    }  
  }

  $script .= <<EOF;
$commandString
date
env
EOF

  `echo '$script' > pbstmp; qsub -V pbstmp; rm pbstmp`;

}

sub launch 
{
  my $self = shift;
  my $semaphore = shift;
  my $logFile = $self -> {STDOUT_LOG_FILE};
  my $logErr = $self -> {STDERR_LOG_FILE};

  my $commandString   = "(" . $self -> {COMMAND} . ") 2>&1 >> $logFile; chmod 664 $logFile; /usr/bin/touch $semaphore;";

  my $mode = $self -> {QUEUE};

  setModuleEnvironment($self->{COMPILER});
  $ENV{MODELERC}=$self->{MODELRC};
  if ($self -> {QUEUE} eq LOCAL) 
  {
    `($commandString) &`; # run in background
  }
  else 
  {
    $self -> runInBatch("$commandString", $self -> {QUEUE});
  }
}

sub setModuleEnvironment 
{
    my $compiler = shift;

    print " setModuleEnvironment: COMPILER=$compiler\n";
    require $ENV{MODELROOT}."/exec/testing/perlreq";

    if ($self->{BRANCH} =~ m/AR5/) 
    {
      if ($compiler eq "intel") 
      {
        module (load, "comp/intel-10.1.017",  "mpi/impi-3.2.2.006");
      } 
      elsif ($compiler eq "gfortran") 
      {
        module (load, "other/comp/gcc-4.5", "other/mpi/openmpi/1.4.2-gcc-4.5");
      } 
      else 
      {
        # Nothing to do
      }
    } 
    else 
    {
      if ($compiler eq "intel")
      {
        module (load, "comp/intel-11.1.072",  "mpi/impi-3.2.2.006");
      }
      elsif ($compiler eq "gfortran")
      {
        module (load, "other/comp/gcc-4.6", "other/mpi/mvapich2-1.4.1/gcc-4.6");
      }
      else 
      {
        # Nothing to do
      }
    }
}

1;