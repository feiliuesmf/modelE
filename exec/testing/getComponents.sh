#!/bin/bash

deck=$1

declare -a c
c=("COMPONENTS :=  ")
compLine="Components:"

filecontent=( `cat $deck ` )

start=0
for entry in "${filecontent[@]}"; do
  if [ "$entry" == "$compLine" ]; then
    start=1
    continue
  fi
  if [ $start -eq 1 ]; then
    if [ "$entry" == "Component" ]; then
      break
    fi
    if [ "$entry" == "tracers" ]; then
      c=( "${c[@]}" "$entry" )
    fi    
  fi    
done
c=( "${c[@]}" " profiler shared solvers giss_LSM dd2d Ent MPI_Support" )

echo ${c[@]} > components
