BEGIN {
    print "set-resources:"
}
/^    [^ ]/ {
    # This matches the main keys (indented 4 spaces)
    if (key != "") {
        # Print all accumulated subkeys for previous main key
        for (subkey in subkeys) {
            print "    - " key ":" subkey "=" subkeys[subkey]
        }
    }
    key = $1
    delete subkeys
    next
}
/^        / {
    # This matches the subkeys (indented 8 spaces)
    subkey = $1
    subkey_value = $2
    subkeys[subkey] = subkey_value
    next
}
END {
    # Print the last accumulated main key's subkeys
    if (key != "") {
        for (subkey in subkeys) {
            print "    - " key ":" subkey "=" subkeys[subkey]
        }
    }
}
