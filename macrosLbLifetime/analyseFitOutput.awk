#!/usr/bin/awk -f

BEGIN {
    # Check arguments
    if (ARGC < 1 ) {
	print "Argument error!";
	    exit(1);
	}

    # Remove argument from ARGV and add implied "-" (stdin)
#    if (ARGC == 1)
#	ARGV[1] = "-";
#    else
#	ARGV[1] = "";
}

/Working on file no/ {
    printf "%s: ", $5;
    next;
}

/Working on/ {
    printf "%s: ", $3;
}

/Fitresults/ {
    print $3, $5, "-", $20, $22, $24, $16, $18;
}

