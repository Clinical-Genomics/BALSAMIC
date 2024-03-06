#!/bin/bash

# Check if at least 3 arguments are provided
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <genome_version> <output_vcf> <tumor_bam> [normal_bam]"
    exit 1
fi

# Assign variables
genome_version="$1"
output_vcf="$2"
tumor_bam="$3"
normal_bam="$4"
output_vcf_tmp="$output_vcf.tmp"

# Print given arguments
echo "genome_version: $genome_version"
echo "output vcf: $output_vcf"
echo "tumor bam: $tumor_bam"
echo "normal bam: $normal_bam"

# Set chr positions depending on the genome version
if [ "$genome_version" = "hg19" ]; then
    igh_chr="14"
    igh_pos="106032614"
    dux4_chr="4"
    dux4_pos="190988100"
elif [ "$genome_version" = "hg38" ]; then
    igh_chr="14"
    igh_pos="105586437"
    dux4_chr="4"
    dux4_pos="190173000"
else
    echo "Invalid genome version. Accepted values: hg19, hg38. Given: $genome_version"
    exit 1
fi


# Define functions
get_supporting_reads() {
      # Get number of supporting reads for IGH::DUX4 rearrangement in a given BAM file
      local bam="$1"
      if [ "$genome_version" = "hg19" ]; then
          local supporting_reads=$(samtools view -F 1024 -c \
           -e '(rnext == "4" && pnext > 190988100 && pnext < 191007000) || (rnext == "10" && pnext > 135477000 && pnext < 135500000) || (rnext == "GL000228.1" && pnext > 70000 && pnext < 115000) || ([SA] =~ "10,1354[789][0-9]{4}") || ([SA] =~ "4,19(09[8-9][0-9]|100[0-7])[0-9]{3}" || [SA] =~ "GL000228.1,([7-9][0-9]{4}|1[0-1][0-5][0-9]{3})")' \
           $bam 14:106032614-107288051 )
      elif [ "$genome_version" = "hg38" ]; then
          local supporting_reads=$(samtools view -F 1024 -c \
           -e '(rnext == "4" && pnext > 190173000 && pnext < 190176000) || ([SA] =~ "4,19017[345][0-9]{3}")' \
           $bam chr14:105586437-106879844 )
      fi
      echo $supporting_reads
}

# Set information for tumor and normal
supporting_reads_tumor=$(get_supporting_reads $tumor_bam)
samples_header="TUMOR"
samples_field="${supporting_reads_tumor}"
if [ -n "$normal_bam" ]; then
    supporting_reads_normal=$(get_supporting_reads $normal_bam)
    samples_header="NORMAL\tTUMOR"
    samples_field="${supporting_reads_normal}\t${supporting_reads_tumor}"
fi


# If supporting reads are found in the tumor, set filter to PASS. Otherwise add: no_supporting_reads
if [ "$supporting_reads_tumor" -gt 0 ]; then
    vcf_filter="PASS"
else
    vcf_filter="no_supporting_reads"
fi

echo "supporting reads tumor: $supporting_reads_tumor"
echo "supporting reads normal: $supporting_reads_normal"
echo "vcf filter: $vcf_filter"

# Write vcf entry
{
  echo '##fileformat=VCFv4.2'
  echo '##ALT=<ID=BND,Description="Break end">'
  echo '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">'
  echo '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">'
  echo '##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of paired-ends that support the event">'
  echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${samples_header}"
  echo -e "${igh_chr}\t${igh_pos}\tsamtools_igh_dux4\tN\tN[${dux4_chr}:${dux4_pos}[\t.\t${vcf_filter}\tSVTYPE=BND;IMPRECISE;\tDV\t${samples_field}"
} >> $output_vcf_tmp

bgzip $output_vcf_tmp > $output_vcf
tabix -p vcf $output_vcf
rm $output_vcf_tmp
