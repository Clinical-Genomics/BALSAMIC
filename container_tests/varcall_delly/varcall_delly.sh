#!/bin/bash -l
# Test if commands exist

valid_commands=( "delly" "bcftools" "tabix" )

for valid_command in "${valid_commands[@]}"
do
  if ! command -v "${valid_command}" &> /dev/null
  then
    echo "${valid_command} could not be found"
    exit 1
  else
    echo "${valid_command} command is found and valid"
  fi
done
