#! /bin/bash -l

# Program to read vcf files and extract the information only 
# for the specified background regions {$1} and output as table
# This output table is required for plotting the detected variants across all samples. For eg: This plot could show the variants trends across samples with different dilutions.

# Usage: bash bcf2AF.sh <background_regions> <input_TNscope_VCF> <input_sample_file_name> <output_table_name> 

source activate D_UMI_APJ

bg_fl=$1
ip_fl=$2
s_nm=$3
op_fl=$4

bcftools query --regions-file $bg_fl -f "%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[%AF\t%AD{0}\t%AD{1}]\n" $ip_fl | awk -v file=$s_nm '{print $1":"$2"_"$3"->"$4"\t"$8/($7+$8)"\t"file}' > $op_fl

