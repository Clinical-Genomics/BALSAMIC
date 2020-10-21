#!/bin/sh
# Bash script adapted from Don Freed (Sentieon). 
# This script allows filtering of consensus reads based on the min_reads settings.

if [ "$#" -ne 3 ] ; then
  echo "Filter duplex UMI consensus reads by group size" >&2
  echo "Usage: $0 INPUT_BAM OUTPUT_BAM MIN_READS" >&2
  exit 1
fi

input_bam=$1
output_bam=$2
min_reads=$3

# min_reads has three values, e.g. 2,1,1
# The first value is applies to the total group size,
# the second value to one single-strand consensus, 
# and the last value to the other single-strand consensus.

if [ -z "$SENTIEON_INSTALL_DIR" ] ; then
    echo "ERROR: Missing SENTIEON_INSTALL_DIR as environment variable in your terminal" >&2
    exit 1
fi

sentieon=$SENTIEON_INSTALL_DIR'/bin/sentieon'

samtools view -h $input_bam | \
awk -v MinR=$min_reads -v OFS="\t" ' 
function min(b) {
   return b[1]>b[2]?b[2]:b[1]
}
function max(b) {
   return b[1]>b[2]?b[1]:b[2]
}
function sum(b) {
   return b[1]+b[2]
}
BEGIN {split(MinR,tmp,","); 
       mr1=tmp[1]; 
       mr2=tmp[2]<tmp[3]?tmp[2]:tmp[3]; 
       mr3=tmp[2]<tmp[3]?tmp[3]:tmp[2];
       reads=0; reads2=0; 
       print "["strftime("Time = %m/%d/%Y %H:%M:%S", systime())"]  Executing FilterDuplexUMIconsensus."> "/dev/stderr"
}
{ if ($0~/^@/) {print;} 
  else {
      for(i=NF;i>=12;i--){ if($i~/^XZ:Z:/) {split($i,a,":");split(a[3],b,","); break;}}
      if (sum(b)>=mr1 && min(b)>=mr2 && max(b)>=mr3 ) {reads2++;print;}
      reads++;
      if (reads % 2000000 == 0) {
          print "["strftime("Time = %m/%d/%Y %H:%M:%S", systime())"]  Processed "reads" reads.  Last read position: "$3":"$4> "/dev/stderr"
      }
  }
}
END {print "["strftime("Time = %m/%d/%Y %H:%M:%S", systime())"]  Done. Processed "reads" reads.  "reads2" reads passed the filter."> "/dev/stderr" } ' | \

samtools view -bh - > $output_bam

$sentieon util index $output_bam
