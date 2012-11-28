#!/usr/local/bin/bash

RESET=$'\033[0m'
BOLD=$'\033[1m'
if [ ! "$2" = "--help" -a ! "$2" = "-h" ]
then
	echo "Usages: $1 [--config myConfig.cfg] or $1 --help for help"; exit 0 ; 
else
echo
echo -e "\t${BOLD}$1${RESET} executes ModelE regression testing."
cat <<pleHHelp1 
During setup, it interprets a number of environmental variables and/or reads a possible configuration
file (see Examples below). The variables will override settings contained in the configuration file.

pleHHelp1
echo -e ${BOLD}Environmental Variables:${RESET}
cat <<pleHHelp2
VARIABLE	EFFECT
--------	------
MODELROOT	Changes to a non-standard root context. Set to \$NOBACKUP/devel by default.
REGWORK		Where to place working files from model run. REGSCRATCH will be appended with /regression_scratch.
		REGRESULTS will be appended with /regression_results.
MODELEBASELINE	Where baseline results are stored for comparison with extended run. 
		Default is $NOBACKUP/modelE_baseline
GCMSEARCHPATH	Data path. Default currently /discover/nobackup/projects/giss/prod_input_files.
	
pleHHelp2
echo -e ${BOLD}Examples:${RESET}
cat <<pleHHelp3
	$1 --config discover.cfg
	Above command would use the preconfigured sectioned config file discover.cfg, which sets all options
	
	$1 
	Above command implies that all needed environment variables are set properly, or that none are set,
	accepting normal default values.
pleHHelp3
fi
