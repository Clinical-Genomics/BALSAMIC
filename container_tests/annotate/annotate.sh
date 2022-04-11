#!/bin/bash
# Test if commands exist

valid_commands=( "bcftools" "vcfanno" "genmod" "vep" "vep_install" )

for tool in ${valid_commands[@]}
do
    tool_msg=$($tool --help 2>&1 | grep 'No')
    if [[ -z $tool_msg ]]; then 
        echo "$tool command exists and works";
    else 
        echo "$tool command doesnt work";
        exit 1;
    fi
done

