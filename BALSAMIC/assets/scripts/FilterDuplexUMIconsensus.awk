#! /bin/awk -f

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
END {print "["strftime("Time = %m/%d/%Y %H:%M:%S", systime())"]  Done. Processed "reads" reads.  "reads2" reads passed the filter."> "/dev/stderr" } 
