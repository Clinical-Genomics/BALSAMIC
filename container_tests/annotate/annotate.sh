#!/bin/bash
# Test if commands exist

valid_commands=( "bcftools" "vcfanno" "genmod" "vep" "vep_install")

for tool in "${valid_commands[@]}"
do
    if ! command -v "${tool}"  &> /dev/null; then
        echo "$tool command not found in the container"
        exit 1;
    else
        tool_msg=$($tool 2>&1 | grep 'No')
        if [[ -z $tool_msg ]];then
            echo "$tool command found and \"$tool --help\" command executes"
        else
            echo "$tool found in the container but \"$tool --help\" command not executing"
            exit 1;
        fi
    fi
done
