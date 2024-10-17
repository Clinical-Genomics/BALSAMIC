#!/usr/bin/awk -f

BEGIN { OFS = "\t" }
/^@/ { print; next }
{
    flag = $2

    # If the mate is unmapped, remove the MC tag if it's "*"
    if (and(flag, 8) != 0) {
        for (i = 12; i <= NF; i++) {
            if ($i ~ /^MC:Z:\*$/) {
                $i = ""
            }
        }
    }

    # Check if any of the specific bitwise flags are set (2, 8, 32, 64, 128)
    if (and(flag, 2 + 8 + 32 + 64 + 128) != 0) {
        # Add mate unmapped flag if the mate is unmapped
        if ($7 == "*" || and(flag, 8) != 0) {
            flag = or(flag, 8)
        }

        # Ensure the read paired flag (1) is set if any of these are present
        flag = or(flag, 1)
    }

    # Set the modified flag
    $2 = flag

    print
}
