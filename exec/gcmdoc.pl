#!/usr/bin/perl
#This script parses fortran files and creates documentation
#in HTML format. It utilizes information provided with @.* tags.
#Author: I. Aleinov
#Version:  1.0
#Usage:  gcmdoc.pl [-D output_dir] file1.f [file2.f ...]
#           output_dir - a directory were to write HTML files

use Getopt::Long;  #module to process command line options

#some global definitions:

$output_dir = "./";

#some useful patterns

$some_decl = "integer|real|double\\s+precision|character|logical|complex";
$some_decl = "(?:$some_decl)(?:\\s*\\*\\s*(?:\\d+|\\(\\d+\\)))?";

print "$some_decl\n";

#an example
#GetOptions("s", "e=s", "f=s", "I=s@", "m=s", "c", "p", "g", "h", "o=s", "a=s")
#                       || die "problem in GetOptions";

GetOptions("D=s") || die "problem in GetOptions";

if( $opt_D ) {
    $output_dir = $opt_D;
    $output_dir =~ s/\/?\s*$//;
}

print "will write output to: $output_dir\n";

#init global hashes:
%db_vars = {};       #var name: module:sub:var
%db_subs = {};       #sub name: module:sub
%db_modules = {};
%db_files = {};


while( $current_file = shift ) {
    open(SRCFILE, $current_file) or die "can't open $current_file\n";
    parse_file();
    close SRCFILE; #just in case want to reset line counter ($.)
}

#print " PRINTING INFO \n";

#while ( ($name,$ref) = each %db_subs ) {
#    print "subroutine: $name\n";
#    print "Comment: $$ref{sum}\n";
#    print "Author: $$ref{auth}\n";
#    print "Version: $$ref{ver}\n";
#}

#print " PRINTING VARIABLES \n";

#foreach $name ( sort keys %db_vars ) {
#    print "$name\n";
#    print "    declaration: $db_vars{$name}{decl}\n";
#    print "    value      : $db_vars{$name}{value}\n";
#    print "    description: $db_vars{$name}{sum}\n"; 
#}

# printing output to html files

print  "printing index\n";
htm_start("index.html","GISS GCM Documentation Index");
print HTM '<H5><Center>GISS GCM Documentation Index</Center></H5>'."\n";
print HTM "<HR>\n";

htm_link("Source Files","files.html"); print HTM "<BR>\n";
htm_link("Fortran Modules","modules.html"); print HTM "<BR>\n";
htm_link("Database Variables","vars_db.html"); print HTM "<BR>\n";
htm_link("Namelist Variables","vars_nl.html"); print HTM "<BR>\n";
htm_link("ALL Variables ( LONG! )","vars_all.html"); print HTM "<BR>\n";

htm_end();

htm_start("files.html","GCM Source Files");
htm_link("Index","index.html");
print HTM '<H3><Center>FILES</Center></font></H3>'."\n";
print HTM "<dl>\n";
foreach $name ( sort keys %db_files ) {
    print HTM "<dt><b>";
    htm_link($name,"$name.html");
    print HTM "</b><BR>\n";
    #print "printing $name\n";
    #print HTM "$name\n";
    #print HTM "<UL>\n";
    #print HTM "Summary: $db_files{$name}{sum}<BR>\n";
    print HTM "<dd>Modules: ";
    #while ( $name_mod = pop @{$db_files{$name}{modules}} ) {
    foreach $name_mod ( @{$db_files{$name}{modules}} ) {
	#print HTM " $name_mod";
	print HTM " ";
	htm_link("$name_mod","$name_mod.html");
    }
    print HTM "<BR>\n";
    #print HTM "Summary: $db_files{$name}{sum}<BR>\n";
    htm_text("Summary: $db_files{$name}{sum}");
}
print HTM "</dl>\n";
htm_end();

htm_start("modules.html","GCM Fortran Modules");
htm_link("Index","index.html");
print HTM '<H3><Center>MODULES</Center></font></H3>'."\n";
print HTM "<dl>\n";
foreach $name ( sort keys %db_modules ) {
    print HTM "<dt><b>";
    htm_link($name,"$name.html");
    print HTM "</b><BR>\n";
    #print HTM "<dd>$db_modules{$name}{sum}<BR>\n";
    print HTM "<dd>";
    htm_text("$db_modules{$name}{sum}");
}
print HTM "</dl>\n";
htm_end();

print  "printing variables\n";

htm_start("vars_db.html","GCM ALL Variables");
htm_link("Index","index.html");
print HTM '<H3><Center>List of Database Variables</Center></font></H3>'."\n";
@var_list = ();
foreach $name ( keys %db_vars ) {
    if ( $db_vars{$name}{tag} =~ /\bdbparam\b/i ) {
	push @var_list, $name;
    }
}
htm_prt_vars( "sl", @var_list );
htm_end();

htm_start("vars_nl.html","GCM ALL Variables");
htm_link("Index","index.html");
print HTM '<H3><Center>List of Namelist Variables</Center></font></H3>'."\n";
@var_list = ();
foreach $name ( keys %db_vars ) {
    if ( $db_vars{$name}{tag} =~ /\bnlparam\b/i ) {
	push @var_list, $name;
    }
}
htm_prt_vars( "sl", @var_list );
htm_end();

htm_start("vars_all.html","GCM ALL Variables");
htm_link("Index","index.html");
print HTM '<H3><Center>List of All Variables</Center></font></H3>'."\n";
htm_prt_vars( "sl", keys %db_vars );
htm_end();

print  "printing files\n";

foreach $name ( keys %db_files ) {
    htm_start("$name.html","$name");
    htm_link("Index","index.html");
    print HTM "<H3><Center>$name</Center></font></H3>\n";
    #$str = $db_files{$name}{sum};
    #$str =~ s/>/&gt;/g; $str =~ s/</&lt;/g; 
    #print HTM "Summary: $str<BR>\n";
    htm_text("Summary: $db_files{$name}{sum}");
    print HTM "Author : $db_files{$name}{auth}<BR>\n";
    print HTM "Version: $db_files{$name}{ver}<BR>\n";
    print HTM "<HR width=10%>\n";
    print HTM "Modules: \n";
    print HTM "<dl>\n";
    foreach $mod_name ( @{$db_files{$name}{modules}} ) {
	print HTM "<dt>";
	htm_link($mod_name,"$mod_name.html");
	print HTM "<br>";
	#print HTM "<dd>$db_modules{$mod_name}{sum}<BR>\n";
	print HTM "<dd>";
	htm_text("$db_modules{$mod_name}{sum}");
    }
    print HTM "</dl>\n";
    print HTM "<HR width=10%>\n";
    print HTM "Global Subroutines: \n";
    print HTM "<dl>\n";
    foreach $sub_name ( sort @{$db_files{$name}{subs}} ) {
	print HTM "<dt>";
	htm_link($sub_name,":$sub_name.html");
	print HTM "<br>";
	$sub_name = ":$sub_name";
	#print HTM "<dd>$db_subs{$sub_name}{sum}<BR>\n";
	print HTM "<dd>";
	htm_text("$db_subs{$sub_name}{sum}");
    }
     print HTM "</dl>\n";
    htm_end();
}

print "printing modules\n";

foreach $name ( keys %db_modules ) {
    htm_start("$name.html","$name");
    htm_link("Index","index.html");
    print HTM "<table width=\"100%\">\n";
    print HTM "<tr align=center>";
    print HTM "<td><H3>$name</font></H3><td>\n";
    print HTM "<td><H5>File: ";
    htm_link(" $db_modules{$name}{file}","$db_modules{$name}{file}.html");
    print HTM "</H5></td>";
    print HTM "</tr>";
    print HTM "</table>\n";
    #$str = $db_modules{$name}{sum};
    #$str =~ s/>/&gt;/g; $str =~ s/</&lt;/g; 
    #print HTM "Summary: $str<BR>\n";
    htm_text("Summary: $db_modules{$name}{sum}");
    print HTM "Author : $db_modules{$name}{auth}<BR>\n";
    print HTM "Version: $db_modules{$name}{ver}<BR>\n";
    #print HTM "File:";
    #htm_link(" $db_modules{$name}{file}","$db_modules{$name}{file}.html");
    print HTM "<HR width=10%>\n";
    print HTM "Subroutines: \n";
    print HTM "<dl>\n";
    #while (  $sub_name = pop @{$db_modules{$name}{subs}} ) {
    foreach $sub_name ( sort @{$db_modules{$name}{subs}} ) {
	print HTM "<dt>";
	htm_link($sub_name,"$name:$sub_name.html");
	print HTM "<br>";
	$sub_name = "$name:$sub_name";
	#print HTM "<dd>$db_subs{$sub_name}{sum}<BR>\n";
	print HTM "<dd>";
	htm_text("$db_subs{$sub_name}{sum}");
    }
    print HTM "</dl>\n";
    print HTM "<HR width=10%>\n";
    print HTM "Global Variables: \n";
    print HTM "<dl>\n";
    foreach $var ( @{$db_subs{"$name:"}{vars}} ) {
	#$var_name = "$name::$var";
	$var_name = $name.":".":".$var;
	#print "VAR::: $var_name\n";
	if ( $db_vars{$var_name}{decl} =~ /parameter/i ) { $color = "#008800" }
	else { $color = "#880000" }
	print HTM "<dt><font color=$color><B>$var</B></font>";
	print HTM " : <code>$db_vars{$var_name}{decl}</code><BR>\n";
	#print HTM "<dd>$db_vars{$var_name}{sum}<BR>\n";
	print HTM "<dd>";
	if ( $db_vars{$var_name}{sum} ) {
	    #print HTM "$db_vars{$var_name}{sum}<BR>\n";
	    htm_text("$db_vars{$var_name}{sum}");
	}
	if ( $db_vars{$var_name}{value} ) {
	  print HTM "Initial Value <code> = $db_vars{$var_name}{value}</code>";
	}
	if ( ! ($db_vars{$var_name}{sum} || $db_vars{$var_name}{value}) ) {
	    print HTM "<BR>\n";
	}
	#print HTM "      <code>$db_vars{$var_name}{decl}</code><BR>\n";
    } 
    print HTM "</dl>\n";
    htm_end();
}


print "printing sub's\n";
foreach $name ( keys %db_subs ) {
    htm_start("$name.html","$name");
    htm_link("Index","index.html");
    #if ( ! $name =~ /([^:]*):([^:]*)/ ) { 
	#print "wrong sub name: $name\n"; next; };
    $sub_name = $name; $sub_name =~ s/^.*://;
    $mod_name = $name; $mod_name =~ s/:.*$//;
    print HTM "<table width=\"100%\">\n";
    print HTM "<tr align=center>";
    print HTM "<td><H3>$sub_name</font></H3><td>\n";
    print HTM "<td><H5>Module: ";
    htm_link("$mod_name","$mod_name.html");
    print HTM "</H5></td>";
    print HTM "<td><H5>File: ";
    htm_link("$db_subs{$name}{file}","$db_subs{$name}{file}.html");
    print HTM "</H5></td>";
    print HTM "</tr>";
    print HTM "</table>\n";


    #print HTM "<H2><Center>$name</Center></font></H2>\n";
    #$str = $db_subs{$name}{sum};
    #$str =~ s/>/&gt;/g; $str =~ s/</&lt;/g; 
    #print HTM "Summary: $str<BR>\n";
    htm_text("Summary: $db_subs{$name}{sum}");
    print HTM "Author : $db_subs{$name}{auth}<BR>\n";
    print HTM "Version: $db_subs{$name}{ver}<BR>\n";
    print HTM "<HR width=10%>\n";
    print HTM "Variables: \n";
    print HTM "<dl>\n";
    #foreach $var_name ( keys %db_vars ) {
	#if ( $var_name !~ /^$name:(\w+)/ ) { next; }
    #while ( $var = pop @{$db_subs{$name}{vars}} ) {
    foreach $var ( @{$db_subs{$name}{vars}} ) {
	$var_name = "$name:$var";
	#print "VAR::: $var_name\n";
	if ( $db_vars{$var_name}{decl} =~ /parameter/i ) { $color = "#008800" }
	else { $color = "#880000" }
	print HTM "<dt><font color=$color><B>$var</B></font>";
	print HTM " : <code>$db_vars{$var_name}{decl}</code><BR>\n";
	print HTM "<dd>";
	if ( $db_vars{$var_name}{sum} ) {
	    #print HTM "$db_vars{$var_name}{sum}<BR>\n";
	    htm_text("$db_vars{$var_name}{sum}");
	}
	if ( $db_vars{$var_name}{value} ) {
	    print HTM 
	      "Initial Value <code> = $db_vars{$var_name}{value}</code><BR>\n";
	}
	if ( ! ($db_vars{$var_name}{sum} || $db_vars{$var_name}{value}) ) {
	    print HTM "<BR>\n";
	}
	#print HTM "      <code>$db_vars{$var_name}{decl}</code><BR>\n";
    }
    print HTM "</dl>\n";
    htm_end();
}



print "done\n";
print "to view the documentation do: \"netscape $output_dir/index.html\"\n";


#### this is the end of the main program ######

sub by_name_sub_mod { 
    @aa = split /:/,$a; 
    @bb = split /:/,$b; 
       lc($aa[2]) cmp lc($bb[2]) ||  
	   lc($aa[1]) cmp lc($bb[1]) || 
	       lc($aa[0]) cmp lc($bb[0])  ;
}

sub htm_prt_vars { # print list of variables
    my $opts = shift;
    my $var_name;
    my $file;
    my $print_location = $opts =~ /l/i;
    my $do_sort        = $opts =~ /s/i;
    my @var_list;
    if ( $do_sort ) {
	@var_list = sort by_name_sub_mod @_;
    } else {
	@var_list = @_;
    }
    print HTM "<dl>\n";
    foreach $var_name ( @var_list ) {
	$var_name =~ /(\w*):(\w*):(\w+)/ || next;
	my $mod = $1;
	my $subr = $2;
	my $var = $3;
	#print "PRT_VARS: $var_name, $mod, $subr, $var\n";

	if ( $db_vars{$var_name}{decl} =~ /parameter/i ) { $color = "#008800" }
	else { $color = "#880000" }
	print HTM "<dt><font color=$color><B>$var</B></font>";
	print HTM " : <code>$db_vars{$var_name}{decl}</code><BR>\n";
	print HTM "<dd>";
	if ( $print_location ) {
	    if ( $subr ) {
		#print HTM "Subroutine: $subr, ";
		print HTM "Subroutine: ";
		htm_link("$subr","$mod:$subr.html");
		$file = $db_subs{"$mod:$subr"}{file};
	    } else {
		print HTM "Global variable ";
		$file = $db_modules{"$mod"}{file};
	    }
	    #print HTM "Module: $mod, File: $file<BR>\n";
	    print HTM ". Module: ";
	    if ( $mod ) { 
		htm_link("$mod","$mod.html");
	    } else {
		print HTM "NONE";
	    }
	    print HTM ". File: ";
	    htm_link("$file","$file.html");
	    print HTM "<BR>\n";
	}
	if ( $db_vars{$var_name}{sum} ) {
	    htm_text("$db_vars{$var_name}{sum}");
	}
	if ( $db_vars{$var_name}{value} ) {
	    print HTM 
	      "Initial Value <code> = $db_vars{$var_name}{value}</code><BR>\n";
	}
	if ( ! ($db_vars{$var_name}{sum} || $db_vars{$var_name}{value}) ) {
	    print HTM "<BR>\n";
	}
	
    }
print HTM "</dl>\n";

}

sub htm_start {
    my $filename = shift;
    my $title = shift;
    open ( HTM, ">$output_dir/$filename" ) || die "can't open $filename";
    print HTM "<html>\n";
    print HTM "<head>\n";
    print HTM "<TITLE>$title</TITLE>\n";
    print HTM "</head>\n";
    print HTM '<body BGCOLOR="#FFFFFF"TEXT="#000000">'."\n";
}

sub htm_end {
    print HTM "</body>\n";
    print HTM "</html>\n";
    close HTM;
}

sub htm_link {
    my $name = shift;
    my $link = shift;
    $link =~ s/:/%3A/g;
    print HTM "<a href=\"$link\">$name</a>";
}

sub htm_text {
    my $str = shift;
    $str =~ s/>/&gt;/g; $str =~ s/</&lt;/g; $str =~ s/&/&amp;/g;
    print HTM "$str<BR>\n";
}

#@sum   UNITNAME Brief summary/description of the current program unit or file
#@auth    Author  
#@ver     Version (a number plus any relevant comments)
#@calls   List of routines called in this program unit
#@cont    List of program units contained within a file or module + other info.
#@fun     FUNCNAME Denotes a function
#@param   PARAMNAME Denotes a parameter (in the FORTRAN sense)
#@var     VARNAME Denotes a variable (in the FORTRAN sense)
#@dbparam VARNAME Denotes a database parameter (in the modelE sense)
#@nlparam VARNAME Denotes a NAMELIST parameter (in the modelE sense)
#@+       Continuation line for @sum/@calls/@cont

sub parse_file {
    print "parsing $current_file\n";
    # resetting globals
    $current_module = "";
    $current_sub = "";
    while( <SRCFILE> ) {
	chop;
	#strip regular comments
	if ( /^C/i || /^![^@]/ ) { next; }
	#parse !@... info here
#	if( /^!\@sum/i ) { #subroutine summary
#	    $doc_tag = "sum";
#	    s/^!\@sum\s+(\w+)\s*//i;
#	    $sum_name = "$current_module:$1";
#	    $current_sub = $1;
#	    $db_subs{$sum_name}{sum} = "$_\n";
#	    next;
#	}


	if( /^!\@(var|param|dbparam|nlparam)\b/i ) { #variable summary
	    $doc_tag = "var";
	    s/^!\@(\w+)\s+((\w+)(\s*,\s*(\w+))*)\s*//i;
	    $tag = $1;
	    #@var_list = split '[ ,]', $1;
	    @var_list = split /\s*,\s*/, $2;
	    while( $name = pop @var_list ) {
		$name =~ tr/A-Z/a-z/;
		$var_name = "$current_module:$current_sub:$name";
		#print "var: $var_name\n";
		#print "comment: $_\n";
		$db_vars{$var_name}{sum} = "$_\n";
		$db_vars{$var_name}{tag} = $tag;
		push @{$db_subs{"$current_module:$current_sub"}{vars}}, $name;
	    }
	    next;
	}

	if( /^!\@(\w+)\b/i ) { #generic tag
	    $doc_tag = $1;
	    s/^!\@\w+\s*//i;
	    if ( $current_sub ) {
		$db_subs{"$current_module:$current_sub"}{$doc_tag} = "$_\n";
		    #print "Sub tag:\n";
		    #print "$current_sub  $doc_tag : $_\n";
	    } elsif ( $current_module ) {
		$db_modules{"$current_module"}{$doc_tag} = "$_\n";
		    #print "Module tag:\n";
		    #print "$current_module  $doc_tag : $_\n";
	    } else {
		$db_files{"$current_file"}{$doc_tag} = "$_\n";
		    #print "File tag:\n";
		    #print "$current_file  $doc_tag : $_\n";
	    }
	    next;
	}

	if( /^!\@\+/i ) { #continuation line
	    #$doc_tag = "var";
	    s/^!\@\+\s*//i;
	    if ( $doc_tag =~ "var" ) {
		while( $name = pop @var_list ) {
		    $var_name = "$current_module:$current_sub:$name";
		    $db_vars{$var_name}{sum} .= "$_\n"; }
	    #} elsif ( $doc_tag =~ "sum" ) {
		#$db_subs{$sun_name}{sum} .= "$_\n";
	    } else { #generic
		if ( $current_sub ) {
		    $db_subs{"$current_module:$current_sub"}{$doc_tag}.="$_\n";
		} elsif ( $current_module ) {
		    $db_modules{"$current_module"}{$doc_tag} .= "$_\n";
		    #print "Module tag:\n";
		    #print "$current_module  $doc_tag : $_\n";
		} else {
		    $db_files{"$current_file"}{$doc_tag} .= "$_\n";
		}
	    }
	    next;
	}
	#print;
	if ( /hjdhsdjhfsjkdf/ ) { print "you are kidding!\n"; }

	if( /^.+!\@(var|param|dbparam|nlparam)\b/i ) { #inline variable summary
	    $doc_tag = "var";
	    $str = $_;
	    $str =~ s/.*!\@(\w+)\s+((\w+)(\s*,\s*(\w+))*)\s*//i;
	    $tag = $1;
	    #$str =~ s/\s*!\@var\s+(.*)$//;
	    #$comment = $1;
	    #$str =~ s/......//;
	    #$str =~ /(
		#      (\w+ (\s*=\s*[^,]+)? \s*,\s*)* \w+ (\s*=\s*[^,]+)?
		#      )\s*$/x;
	    @var_list = split /\s*,\s*/, $2;
	    while( $name = pop @var_list ) {
		$name =~ tr/A-Z/a-z/;
		$var_name = "$current_module:$current_sub:$name";
		$db_vars{$var_name}{sum} = "$str\n";
		$db_vars{$var_name}{tag} = $tag;
		push @{$db_subs{"$current_module:$current_sub"}{vars}}, $name;
	    }
	    s/!\@.*$//i;
	}

	# prarse fortran code here
	if ( /!/ ) {  # possible comment at the en of the line
	    /^(
	       (  [^"'!] | (\"[^"]*\") | (\'[^']*\')  )*
              )/x;  # '))) / ; -emacs fix
	    $_ = $1;
	    #print "$_\n";
	}

	if (  /^     \S/ ) { # continuation line
	    s/^     \S/ /;
	    $fstr .= $_;
	    next;
	}

	# parse previous line
	parse_fort_str();

	#first check if beginning/end of block

	if ( /^\s*module\s+(\w+)/i ) { #start module
	    if ( $1 !~ /procedure/i ) {
		$current_module = $1;
		$current_sub = "";
		$db_modules{"$current_module"}{file} = $current_file;
		push @{$db_files{"$current_file"}{modules}}, $current_module;
		#print "MOD: ", ${$db_files{"$current_file"}{modules}}[0],"\n";
	    }
	}
	if ( /^\s*end\s+module\b/i ) { #end module
	    $current_module = "";
	}
	if ( /^\s*(subroutine|(?:$some_decl\s+)?function|program)\s+(\w+)/i ) {
	    $current_sub = $2;
	    $db_subs{"$current_module:$current_sub"}{file} = $current_file;
	    if ( $current_module ) {
		push @{$db_modules{"$current_module"}{subs}}, $current_sub;
	    } else {
		push @{$db_files{"$current_file"}{subs}}, $current_sub;
		#push @{$db_modules{"\@GLOBAL"}{subs}}, $current_sub;
	    }
	}
	if ( /^\s*end\s+(subroutine|function|program)\b/i 
	     || /^\s*end\s*$/i ) {
	    $current_sub = "";
	}
	if ( /^\s*(block\s+data)\s+(\w+)/i ) { # block data hack
	    $current_sub = 'BLOCK_DATA_'.$2;
	}

	# beginning of the new fortran line
	$fstr = $_;
				   
    }
    print "number of lines: $.\n";
}

sub parse_fort_str {
    #print "$fstr\n";

    # variable declaration string
    if ( $fstr =~ /^\s*(integer|real|double|character|logical|complex)/i ) {
	if ( $fstr =~ s/^\s*(\S+.*)::\s*//i ) {
	    $var_type = $1;
	} else {
	    $fstr =~ s/^\s*(\S+)\s*//i;
	    $var_type = $1;
	    if ( $var_type =~ /double/i ) {
		$var_type = "real*8";
		$fstr =~ s/^\s*(\S+)\s*//i;
	    }
	}
	$var_type =~ s/ //g;
	$var_type =~ tr/A-Z/a-z/;
	#print "VAR_TYPE:  $var_type\n";
	# parse single variables
	while ( $fstr =~ s/^\s*,*\s*(\w+)((\([^()]+\))*)// ) {
	    $var = $1;
	    $dim = $2;
	    $var =~ tr/A-Z/a-z/;
	    #print "VAR $var   DIM $dim\n";
	    # look for initial value
	    $var_val = "";
	    # var = abc( ()), ( ())
	    if( $fstr =~ 
		s/^\s*=\s*(
			   [^,"'()]* (   \(  ( [^()] | (\([^()]*\)) )*  \)  )
                          )//x ) {
              #  " ))) {
		$var_val = " $1";
		#print "VAR:  $var   ::: $var_val\n";
	    }
	    # var = "abc", 'abc', 123
	    if( $fstr =~ s/^\s*=\s*(\"[^"]*\"|\'[^']*\'|[^,"'()]+)// ) {#'))) {
		#print "FSTR: #$fstr#\n";
		#if( $fstr =~ s/^\s*=\s*(\".+\")// ) {
		#print "VAL: #$1#\n";
		$var_val = " $1";
		#print "VAR:  $cvar   ::: $var_val\n";
	    }
	    $var_name = "$current_module:$current_sub:$var";
	    $db_vars{$var_name}{decl} = $var_type;
	    if ( $dim ) { $db_vars{$var_name}{decl} .= ",dimension$dim"; }
	    $db_vars{$var_name}{value} = $var_val;
	    if ( ! $db_vars{$var_name}{sum} ) { #wasn't included earlier
		push @{$db_subs{"$current_module:$current_sub"}{vars}}, $var;
	    }
	    #if ( ! $current_sub ) {
		#print "GLOBAL: $current_module $var_name\n";
	    #}
	}
    return;
    }
}

sub old_old {
while ($filename = shift) {
  open(IN, $filename) or die "can't open $filename\n";
  #print "reading $filename\n";
  $str = "";
  $function = "_GLOBAL_";
  $module = "_NONE_";
  $print_header = 1;
  while (<IN>) {
    if ( /^[Cc]/ || /^ *!/ || /^\s*$/ ) {next;}
    if (  /^     \S/ ) {
      #print "cont: $_\n";
      $str .= $_;
      #print;
    }
    else {
      if ( $str =~ /^\s*(module)\s+(\w+)/i ) {
        $module = $2;
        $function = "_GLOBAL_";
        #print "FUCTTION:: $function\n";
        $print_header = 1;
      }
      if ( $str =~ /^\s*(subroutine|function|program)\s+(\w+)/i ) {
        $function = $2;
        #print "FUCTTION:: $function\n";
        $print_header = 1;
      }
      if ( $str =~ /^\s*(block\s+data)\s+(\w+)/i ) {
        $function = 'BLOCK_DATA_'.$2;
        #print "FUCTTION:: $function\n";
        $print_header = 1;
      }
      if ( $str =~ /\b$var\b/i ) { 
        $is_decl = ($str =~ /^\s*(integer|real|double|logical|character|data|dimension)/i);
        $is_used = ($str =~ /^\s*(use|common)/i);
        $is_changed = ($str =~ /\b$var\b\s*(\(.*\))?\s*=|call.*\b$var\b/i);
        if( ($show_decl && $is_decl) || 
            ($show_use && $is_used) ||
            ($show_body && ( (!($is_decl||$is_used)) || $str=~/=> *\b$var\b/i))
          ||($show_changed && $is_changed) ) {
          if ( $print_header ) {
            print "\nFILE:: $filename  MODULE:: $module  SUB:: $function\n\n";
            $print_header = 0;
          }
          $str =~ s/\b($var)\b/&to_bold($1)/egi;
          print "$. :\n";
          print $str;
        }
      }
      $str = $_;
    }
  }
}

} # sub old_old

sub to_bold {
  my $str = shift(@_);
  #print "got $str\n";
  $str =~ s/(.)/$1\x08$1/g;
  #print "made: $str\n";
  return $str;
}

