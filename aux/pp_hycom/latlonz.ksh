#!/usr/bin/ksh
USAGE="$0  arg1 ..."
# SCRIPT: NAME_of_SCRIPT
# AUTHOR: Nick Tausnev, ntausnev@giss.nasa.gov
# DATE:   DATE_of_CREATION 7/9/2010 
# REV:    1.1.A (Valid are A, B, D, T and P)
#               (For Alpha, Beta, Dev, Test and Production)
#
# PLATFORM: (SPECIFY: AIX, HP-UX, Linux, Solaris 
#                      or Not platform dependent)
#
# PURPOSE: Give a clear, and if necessary, long, description of the
#          purpose of the shell script. This will also help you stay
#          focused on the task at hand.
#
# REV LIST:
#        DATE:  DATE_of_CREATION 7/9/2010 
#        BY:	  AUTHOR_of_MODIFICATION   
#        MODIFICATION: Describe what was modified, new features, etc-
#
#
# set -n   # Uncomment to check your syntax, without execution.
#          # NOTE: Do not forget to put the comment back in or
#          #       the shell script will not execute!
# set -x   # Uncomment to debug this shell script (Korn shell only)
#          
##########################################################
########### DEFINE FILES AND VARIABLES HERE ##############
##########################################################


# HARD CODING need change later 
#latlonz_exe=/discover/nobackup/ntausnev/HYCOM_DIAG/aux/pp_hycom/latlonz
latlonz_exe=./latlonz

##########################################################
############### DEFINE FUNCTIONS HERE ####################
##########################################################

function help_use
{
#
# Display help message and quit
#
cat << ENDOFTEXT

Dear $USER, the usage of the script $0 is as follows:
usage:   $0 [-h]  [-o outFile  -t "title" -i file1 [ file2 [file*] ]
example: ( input files can be zip files !!! )
   $0 \\
     -o /discover/nobackup/ntausnev/Z_Deleted_Files/JAN2280_82.zoutEhMay2 \\
     -t "RunId=EhMay2 Mon=JAN YEAR=2280-82" \\
     -i /discover/nobackup/ntausnev/RUNS_ME/EhMay2/00HYC/outEhMay2_2280-2289/JAN228[0-2].outEhMay2.gz

ENDOFTEXT
exit 1
}

##########################################################
################ BEGINNING OF MAIN #######################
##########################################################

if (( $# < 3 )) 
then
    help_use
fi

vflag=off
filename=
title=

while getopts ht:o:i: opt
do
    case "$opt" in
      h)  hflag=on ; help_use ;;
      t)  title="$OPTARG";;
      o)  fileOUT="$OPTARG";;
      i)  fileIN="$OPTARG";;    # takes first file 
      \?)       # unknown flag
      print >&2  "usage: $0 [-h] [-o outFile -t \"title\"  -i file1 [file2 ..] ] "
          exit 1;;
    esac
done
shift `expr $OPTIND - 1`

listArg="$fileIN $@" #  takes ALL input files
print "\nInput list of files for processing (from command line):"
print $listArg | awk '{ for ( i = 1; i <= NF; i++ ) print $i}'

list=$(ls $listArg 2>/dev/null )
print "\nNext files are found on disk :\n$list\n"

# The first check ALL input files are exist and readable
for ifile in $list
do
  if [[ ! -r  $ifile && ! -r ${ifile}.gz ]]
  then
    print "\nInput file=${ifile} or .gz  does not exist or not readable"
    print "exit from script=$0 with exit code=1"
    exit 1
  fi
done
 

# If some files are gzip then gunzip at ___work directory
rm -rf ___work 2> /dev/null
mkdir ___work

#Copy (and gunzip if need) input files at ___work directory
files=''
for ifile in $list
do 
  if [[ $ifile != *.gz ]]
  then 
     files="${files} ${ifile}"
  else
     cp -p ${ifile} ___work/. ; nfile="___work/`basename ${ifile%.gz}`"
     gunzip ${nfile}.gz ; files="${files} ${nfile}"
  fi
done

command="${latlonz_exe} $fileOUT \"$title\"  $files"
print "\nExecution command:\n   $command"
eval $command
rm -rf ___work 2> /dev/null

print "\n\nScript: $0 ended"
exit 0
# End of script

