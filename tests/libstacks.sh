# Stacks Test Suite library.
# 
# Provides test functions for the Stacks Test Suite.  The prominent functions
# (e.g. ok_() and skip_()) are cousins of similar functions found in autoconf's
# libtap.sh.
#
# Functions prepended with '_' are internal functions not for use by test
# engineers.

source $test_path/libtap.sh

# Return the list of files which are missing from output (as compared to the
# expected output).
# Arguments:
#     eout_path - The directory containing the expected output files.
#     out_path  - The directory containing the actual output files.
_compare_output_file_names () {
    eout_path="$1"
    out_path="$2"
    missing_files=""

    # Check arguments
    if [ -z "$eout_path" ] || [ -z "$out_path" ] ; then
        echo "ERROR: $FUNCNAME() takes exactly two arguments."
    fi

    # Compare expected and actual output files
    for F in `/bin/ls $eout_path`; do
        if [ ! -f $out_path/$F ] ; then
          missing_files="$missing_files $F"
        fi
    done
    echo "$missing_files"
}

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

        # First, determine if this is gzipped or not.
        filename=$(basename "$F")
        ext="${filename##*.}"
        if [ "$ext" = 'gz' ] || [ "$ext" = 'gzip' ] ; then
            test_compression=true
            diff_prog='zdiff'
        else
            test_compression=false
            diff_prog='diff'
        fi

        # If file is expected to be compressed, confirm that it actually is
        # compressed.  Use case we're trying to test for is when $out_path/a.gz
        # is not compressed (despite its misleading 'gz' extension).
        if [ "$test_compression" = true ] ; then
            if [[ `gzip -t $out_path/$F 2>&1` ]] ; then
                echo "Strange error. $out_path/$F actually is NOT compressed."
            fi
        fi

        diff_cmd="$diff_prog -Naur $eout_path/$F $out_path/$F"
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

    cmd=`echo $cmd | sed -e "s@%in@$in_path@"`
    cmd=`echo $cmd | sed -e "s@%out@$out_path@"`

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
# Flags:
#       -i path  - Indicate that output will be 'inline' with input, and that 
#                  the input directory is found at 'path'.  In other words, the
#                  stacks tool under tests places its output in the same
#                  directory as its input.  Behind the scenes, the test harness
#                  makes a copy of 'path' to $out_path/ and sets
#                  in_path=$out_path.  You *must* use '%in' in the 'cmd'
#                  argument.  You may use -i '%in' if the input may be found in
#                  the standard location.
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
    local OPTIND
    bad_flag=false
    inline_path=''

    while getopts ":i:" opt; do
        case "${opt}" in
            i)
                inline_path="${OPTARG}"
                ;;
            *)
                bad_flag=true
                ;;
        esac
    done
    shift $((OPTIND-1))

    tap_desc="$1"
    io_path="$2"
    cmd="$3"

    eout_path="$test_data_path/$io_path/eout"
    in_path="$test_data_path/$io_path/in"
    
    # If -i '%in', then fixup inline_path.
    if [ "$inline_path" = '%in' ] ; then 
        inline_path="$in_path"
    fi

    # Check args
    if [ "$bad_flag" = true ] ; then 
        # User passed a bad flag to ok_()
        echo not ok "$count$tap_desc - ERROR: Bad flag sent to $FUNCNAME."
        failed=`expr $failed + 1`
    elif [ "$inline_path" ] && [ ! -d "$inline_path" ] ; then
        # inline_path is not a directory
        echo not ok "$count$tap_desc - ERROR: $inline_path not a directory in $FUNCNAME."
        failed=`expr $failed + 1`
    elif [ -z "$tap_desc" ] || [ -z "$io_path" ] || [ -z "$cmd" ] ; then
        # These crucial variables had better be set.
        echo not ok "$count$tap_desc - ERROR: Bad call of $FUNCNAME."
        failed=`expr $failed + 1`
    else

        # Setup 
        if [ -n "$inline_path" ] ; then 
            # Tool puts output in input directory (i.e. inline).
            # Copy $inline_path contents to the temp outdir and update $in_path.
            cp -a $inline_path $out_path
            in_path="$out_path"
        else
            # Normal operation.
            mkdir $out_path
            in_path="$test_data_path/$io_path/in"
        fi

        if [ -n "$tap_desc" ] ; then
            tap_desc=" - $tap_desc"
        fi

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
            echo "  ## The following command failed with non-zero status:"
            echo "  $cmd"
            echo "$tap_output" | sed 's/^/  /'
        else

            # Confirm file contents the same...
            file_names_missing="$(_compare_output_file_names $eout_path $out_path)"
            diff_output="$(_compare_output_files $eout_path $out_path)"

            if [ -z "$file_names_missing" ] && [ -z "$diff_output" ]; then
                echo ok "$count$tap_desc"
            elif [ -n "$file_names_missing" ] ; then
                failed=`expr $failed + 1`
                echo not ok "$count$tap_desc"
                echo "  ## The following command failed with missing files:"
                echo "  $cmd"
                echo "$tap_output" | sed 's/^/  /'
                echo "  ## Files missing are:"
                echo "$file_names_missing" | sed 's/^/  /'
                echo "  ## Files in output directory are:"
                echo "$(ls -l $out_path)" | sed 's/^/  /'
            else
                failed=`expr $failed + 1`
                echo not ok "$count$tap_desc"
                echo "  ## The following command failed with unexpected output:"
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
