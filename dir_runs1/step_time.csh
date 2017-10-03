#!/bin/csh
# processing the log file using awk file stats.awk

set infile = $1

if ("$infile" == "") then
  set infile = log
endif
if ( ! -f $infile ) then
  echo "ERROR: file $infile not found"
  exit
endif

if (-f stats.awk) then
  grep '^cpu at step' $infile | awk '{print $8}' | awk -f stats.awk
else
  if (-f ../stats.awk) then
    grep '^cpu at step' $infile | awk '{print $8}' | awk -f ../stats.awk
  else
    echo "stats.awk not found and ../stats.awk not found"
  endif
endif
