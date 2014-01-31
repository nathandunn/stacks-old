# Stacks Test Suite library.
# 
# Provides test functions for the Stacks Test Suite.  The prominent functions
# (e.g. ok_() and skip_()) are cousins of similar functions found in autoconf's
# libtap.sh.
#
# Functions prepended with '_' are internal functions not for use by test
# engineers.

source $test_path/libtap.sh

# Return the differences between expected and actual output files.
# Arguments:
#     eout_path - The directory containing the expected output files.
#     out_path  - The directory containing the actual output files.
_compare_output_files () {
    eout_path="$1"
    out_path="$2"

    # Check arguments
    if [ -z "$eout_path" ] || [ -z "$out_path" ] ; then
        echo "ERROR: $FUNCNAME() takes exactly two arguments."
    fi

    # Compare expected and actual output files
    for F in `/bin/ls $eout_path`; do
        diff_cmd="diff -Naur $eout_path/$F $out_path/$F"
        diff_output="$($diff_cmd)"
        if [ -n "$diff_output" ]; then
            echo "$diff_cmd"
            echo "$diff_output"
        fi
    done
}

# Interpolate the command string to replace %out and %in with the correct full
# path directories.
_interpolate_command () {
    cmd="$1"
    in_path="$2"
    out_path="$3"

    cmd=`echo $cmd | sed -e "s@[[:blank:]]%in[[:blank:]]@ $in_path @"`
    cmd=`echo $cmd | sed -e "s@[[:blank:]]%out[[:blank:]]@ $out_path @"`

    echo $cmd
}

# Mimic of libtap's skip().
# Skip the next test.  Takes the reason why the test is skipped.
# Arguments: 
#       tap_desc - A string containing the test description.
#       io_path  - A string containing the full path to the input and expected
#                  output for the command under test. 
#       cmd      - A string containing the command under test.
# Required global variables:
#       $count          - This variable is used by libtap to track the number
#                         of tests run so far.
#       $failed         - This variable is used by libtap to track the number
#                         of tests that have failed so far.
#       $out_path       - This variable contains the directory to which the
#                         command under test will place its output files. 
#       $test_data_path - This variable contains the directory which is the
#                         parent directory of $io_path.
skip_ () {
    tap_desc="$1"
    io_path="$2"
    cmd="$3"

    in_path="$test_data_path/$io_path/in"

    cmd=$(_interpolate_command "$cmd" $in_path $out_path)
    
    echo "ok $count # skip $tap_desc: $cmd"
    count=`expr $count + 1`
}

# Mimic of libtap's ok().
# Arguments: 
#       tap_desc - A string containing the test description.
#       io_path  - A string containing the full path to the input and expected
#                  output for the command under test. 
#       cmd      - A string containing the command under test.
# Required global variables:
#       $count          - This variable is used by libtap to track the number
#                         of tests run so far.
#       $failed         - This variable is used by libtap to track the number
#                         of tests that have failed so far.
#       $out_path       - This variable contains the directory to which the
#                         command under test will place its output files. 
#       $test_data_path - This variable contains the directory which is the
#                         parent directory of $io_path.
ok_ () {
    tap_desc="$1"
    io_path="$2"
    cmd="$3"
    
    in_path="$test_data_path/$io_path/in"
    eout_path="$test_data_path/$io_path/eout"

    # Check args
    if [ -z "$tap_desc" ] || [ -z "$io_path" ] || [ -z "$cmd" ] ; then
        # TODO - replace with 'bail out'
        echo not ok "$count$tap_desc - ERROR: Bad call of $FUNCNAME."
        failed=`expr $failed + 1`
    else
        # Setup
        if [ -n "$tap_desc" ] ; then
            tap_desc=" - $tap_desc"
        fi
        mkdir $out_path

        # Fix up command
        cmd=$(_interpolate_command "$cmd" $in_path $out_path)

        # Run the command under test
        tap_output=`eval $cmd 2>&1`
        tap_status=$?

        # Check output
        if [ 0 != $tap_status ] ; then 
            # Bad return status.
            failed=`expr $failed + 1`
            echo not ok "$count$tap_desc"
            echo "$tap_output" | sed 's/^/  /'
        else
            # Confirm files same...
            diff_output="$(_compare_output_files $eout_path $out_path)"
            if [ -z "$diff_output" ]; then
                echo ok "$count$tap_desc"
            else
                failed=`expr $failed + 1`
                echo not ok "$count$tap_desc"
                echo "  ## The following command failed:"
                echo "  $cmd"
                echo "$tap_output" | sed 's/^/  /'
                echo "  ## Differences between expected and actual output:"
                echo "$diff_output" | sed 's/^/  /'
            fi
        fi
        count=`expr $count + 1`
        rm -rf $out_path
    fi
}
