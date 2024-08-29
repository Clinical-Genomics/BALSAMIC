#!/usr/bin/awk -f

BEGIN {
    ENVIRON["LC_ALL"] = "en_US.UTF-8"
}

# If the line starts with a '#', it's a header, so print it as is
$1 ~ /^#/ {
    print $0;
    next;
}

# Otherwise, send the body lines to an external sort command
{
    print $0 | "/usr/bin/sort -k1,1V -k2,2n"
}