#!/usr/bin/perl

while ($_ = $ARGV[0], /^-/) {
  shift;
  last if /^--$/;
  if (/^-par\b/) { $rfile = shift; next;}
}

print "parfile=$rfile\n";

open(PARFILE, "$rfile") or die "can't open $rfile";

$exist=1;

while (<PARFILE>) {
    if ( /End/ ) { last; }
    s/!.*//;
    push @parameter, /([\w.+_-]+\s*=\s*[\w.\/+_-]+)/g;
}

foreach $_ ( @parameter ) {
  ($name, $dest) = split /\s*=\s*/;
  if ( $name !~ /^(filedir|filesource|filetarget|regridfile|imsource|jmsource|ntilessource|imtarget|jmtarget|ntilestarget)$/ ) {
    print "parameter $name doesn't exist\n";
    $exist=0
    #exit 1;
  } else {
    if ($name =~ "filesource") {
      $filesource=$dest;
    }
    elsif ($name =~ "filetarget") {
      $filetarget=$dest;
    }
    elsif ($name =~ "filedir") {
      $filedir=$dest;
    }
    elsif  ($name =~ "regridfile") {
      $regridfile=$dest;
    }
    elsif  ($name =~ "imsource") {
      $imsource=$dest;
    }
    elsif  ($name =~ "jmsource") {
      $jmsource=$dest;
    }
    elsif  ($name =~ "ntilessource") {
      $ntilessource=$dest;
    }
    elsif  ($name =~ "imtarget") {
      $imtarget=$dest;
    }
    elsif  ($name =~ "jmtarget") {
      $jmtarget=$dest;
    }
    elsif  ($name =~ "ntilestarget") {
      $ntilestarget=$dest;
    }

  }

}

if ( $exist) {
  $files1="$filedir/$filesource";
  print "$files1\n";
  `cp $files1 .`;
  $rfile1="$filedir/$regridfile";
  print "$rfile1\n";
  `cp $rfile1 .`;

  `./ll2cs $filesource $filetarget $regridfile $imsource $jmsource $ntilessource $imtarget $jmtarget $ntilestarget > output`;

  open( FILE, "< output" ) or die "Can't open output file : $!";
  while( <FILE> ) {
    print;
  }
  close FILE;
}
