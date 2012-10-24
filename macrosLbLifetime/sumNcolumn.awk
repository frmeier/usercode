#!/usr/bin/awk -f

# Calculate sum of the Nth column
# Usage sumNcolumn.awk

BEGIN {
	# Check arguments
	if (ARGC < 2 || int(ARGV[1]) == 0) {
		    print "Argument error!";
			    exit(1);
				}

	    # Use argument to set N
	    N = ARGV[1];
		
		# Remove argument from ARGV and add implied "-" (stdin)
		if (ARGC == 2)
			    ARGV[1] = "-";
		    else
				ARGV[1] = "";

# initialise for sum
			s = 0;
			    n = 0;
}

{
        s+=$(N);
	    n++;
}

END {
        print "Sum of column no. " N ": ", s, " (", n, "rows counted)"
}

