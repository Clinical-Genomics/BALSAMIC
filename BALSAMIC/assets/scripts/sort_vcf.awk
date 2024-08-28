#!/usr/bin/awk -f

# If the line starts with a '#', it's a header, so print it as is
$1 ~ /^#/ {
    print $0;
    next;
}

# Otherwise, send the body lines to an external sort command
{
    print $0 | "sort -k1,1V -k2,2n"
}