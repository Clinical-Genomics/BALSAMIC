# Prints column headers from ucsc .sql schema files

BEGIN{ OFS="\t" }

/CREATE TABLE/ { flag=1; next }
/ENGINE/ { flag=0 }
flag && $0!~/KEY/ {
  gsub(/^[ ]+/,"",$0);
  gsub("`","",$0);
  header=( header=="" ? $1 : header OFS $1)
}END{
  print header
} 
