#!/usr/bin/python 

# This script accepts as an argument the name of a file containing
# fixed format Fortran source code and attempts to create a similarly
# named file containing free format Fortran source code.  By
# convention, the fixed format file has th e '.f' suffix and the free
# format file has '.F90'.

# There is are two situations that this script willl process incorrectly.
#
# Case 1:  A fortran token contains an internal space character.  In fixed 
#          format, spaces are not signficant and the following line is legal:
#               if (x .ne . y) then
#                        ^ - interior space
#          In free format it must be expressed as:
#               if (x .ne. y) then
# Case 2:  Continuations that are sandwiched by #ifdef constructs can be wrong.


import sys
import os
import re
import tempfile

emptyRegexp = re.compile('^(\s*)\n')
commentRegexp= re.compile('^([!cC#]|\s+!)',re.IGNORECASE)
continueRegexp = re.compile('^     \S')
trailingCommentRegExp = re.compile('!')

def markWillContinue(line):
    if trailingCommentRegExp.search(line):
        statement = line.replace('!','& !')
    else:
        statement = line.rstrip() + ' &\n'
    return statement

def markContinuesPrevious(line):
    return '     &' + line[6:]

def convertToFreeFormat(lines):
    newLines = []
    lastStatement = -1

    j = 0;
    for line in lines:
        if (commentRegexp.match(line) or emptyRegexp.match(line)):
            newline = line
            if newline[0] == 'c' or newline[0]== 'C':
                newline = '!' + newline[1:]
            newLines.append(newline)
        else:
            if (continueRegexp.match(line)) :
                previousStatement = newLines[lastStatement]
                newLines[lastStatement] = markWillContinue(previousStatement)
                newLine = markContinuesPrevious(line)
            else:
                newLine = line
            newLines.append(newLine)
            lastStatement = j
        j = j + 1
    return newLines


inFileName = sys.argv[1]
outFileName = inFileName.replace('.f','.F90')

inFile = open(inFileName,"r")
outFile = open(outFileName,"w")

for line in convertToFreeFormat(inFile.readlines()):
    outFile.write(line)

inFile.close()
outFile.close()


