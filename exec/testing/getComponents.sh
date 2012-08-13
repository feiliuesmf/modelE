#!/bin/bash

deck=$1

declare -a c
c=("COMPONENTS := ")
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
    c=( "${c[@]}" "$entry" )
  fi    
done

echo ${c[@]} > components
