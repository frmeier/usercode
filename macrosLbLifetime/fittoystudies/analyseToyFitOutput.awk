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

/Thisplot/ {
    printf "%s ", $2;
}

/Mean/ {
    printf "%s %s ", $3, $4;
}

/Sigma/ {
    printf "%s %s\n", $3, $4;
}

