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

BEGIN { tau=0; tauE=0; nsig=0; npr=0; nnpr=0;}

/ tau / { tau=$3; tauE=$4 }
/ n_signal / { nsig = $3 }
/ n_prompt / { npr = $3 }
/ n_nonprompt / { nnpr = $3 }

/Working on/ { print tau, tauE, nsig, npr, nnpr}

END { print tau, tauE, nsig, npr, nnpr}

