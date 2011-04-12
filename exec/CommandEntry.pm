#!/usr/bin/perl
package CommandEntry;
use Env;

sub new {
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

  # override defaults with arguments
  while ( my ($key, $value) = each(%$args) ) {
    $self -> {$key} = $value;
  }

  bless ($self, $class);
  return $self;
}

sub dependency {
  my $self       = shift;

  return ($self -> {DEPENDENCY});
}

sub status {
  my $self = shift;
  return ($self -> {STATUS});
}

sub isComplete {
  my $self       = shift;

  return ( $self -> {STATUS} eq COMPLETE );
}

sub isRunning {
  my $self       = shift;

  return ( $self -> {STATUS} eq RUNNING );
}

sub isPending {
  my $self       = shift;

  return ( $self -> {STATUS} eq PENDING );
}

sub isReady {
  my $self       = shift;
  
  my $pending = ($self -> isPending());
  my $dependency = ($self -> dependency());

  if ($dependency eq "") { return $pending }
  else {return ($pending and ($self -> dependency() -> isComplete()) )};
}

sub setStatus {
  my $self = shift;
  my $newStatus = shift;
  $self -> {STATUS} = $newStatus;
}

sub setDependency {
  my $self = shift;
  my $newDependency = shift;
  $self -> {DEPENDENCY} = $newDependency;
}

sub runInBatch {
  my $self = shift;
  my $commandString = shift;
  my $queue = shift;
  my $queueString = " ";
  if ($queue) {$queueString = "-q $queue\n";};

#  my $NODE_SIZE = 8;
  my $NODE_SIZE = 16;

  my $nCPUS = $self -> {NUM_PROCS};
  my $nodes = 1 + int(($nCPUS-1)/$NODE_SIZE);
  my $ncpus = $nCPUS;
  $ncpus = $NODE_SIZE if ($ncpus > $NODE_SIZE);

  $walltime = "3:00:00\n";

  print " runInBatch: COMPILER=$self->{COMPILER}\n";
  my $script = <<EOF;
#!/bin/bash
#PBS -l select=$nodes:ncpus=8:proc=neha
#PBS -l walltime=$walltime
#PBS -W group_list=a940a
#PBS -N regression
#PBS -j oe
#PBS $queueString

cd \$PBS_O_WORKDIR


. /usr/share/modules/init/bash
module purge
EOF

    if ($self->{COMPILER} eq intel) {

  $script .= <<EOF;
module load comp/intel-11.1.072  mpi/impi-3.2.2.006
EOF

} else {

  $script .= <<EOF;
module load other/comp/gcc-4.5
export PATH=/usr/local/other/openMpi/gcc-4.5/bin:\$PATH
export LD_LIBRARY_PATH=/usr/local/other/openMpi/gcc-4.5/lib:\$LD_LIBRARY_PATH
EOF

}

$script .= <<EOF;
$commandString
date
env
EOF

`echo '$script' > pbstmp; qsub -V pbstmp; rm pbstmp`;

}

sub launch {
  my $self = shift;
  my $semaphore = shift;
  my $logFile = $self -> {STDOUT_LOG_FILE};
  my $logErr = $self -> {STDERR_LOG_FILE};

  my $commandString   = "(" . $self -> {COMMAND} . ") 2>&1 >> $logFile; /usr/bin/touch $semaphore;";

  my $mode = $self -> {QUEUE};

  print " launch: COMPILER=$self->{COMPILER}\n";
  setModuleEnvironment($self->{COMPILER});
  $ENV{MODELERC}=$self->{MODELRC};
  if ($self -> {QUEUE} eq LOCAL) {
    `($commandString) &`; # run in background
  }
  else {
    $self -> runInBatch("$commandString", $self -> {QUEUE});
  }
}

sub setModuleEnvironment {
    my $compiler = shift;
    print " setModuleEnvironment: COMPILER=$compiler\n";
    require "$ENV{MODULESHOME}/init/perl";  
    module (purge);
    if ($compiler eq intel) {
	module (load, "comp/intel-11.1.072",  "mpi/impi-3.2.2.006");
    } elsif ($compiler eq gfortran) {
	module (load, "other/comp/gcc-4.5");
	$ENV{LD_LIBRARY_PATH}="/usr/local/other/openMpi/gcc-4.5/lib:".$ENV{LD_LIBRARY_PATH};
	$ENV{PATH}="/usr/local/other/openMpi/gcc-4.5/bin:".$ENV{PATH};
    } else {
    }
}

1;
