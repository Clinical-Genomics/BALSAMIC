#!/bin/bash -l
# Test if commands exist

valid_commands=( "bcftools" "samtools" "tabix" "vardict" "vardict-java" )

for valid_command in "${valid_commands[@]}"
do
  if ! command -v "${valid_command}" &> /dev/null
  then
    echo "${valid_command} could not be found"
    exit 1
  fi
done
