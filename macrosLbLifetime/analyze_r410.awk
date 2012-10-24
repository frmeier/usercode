BEGIN {
    ls = 0;
    ev = 0;
    prob = 0;
    mlb = 0;
    mcmatch = 0;
}

{
    if ($4==ls && $6==ev)
    {
	if (prob>$8)
	    print ls, ev, prob, mlb, mcmatch," ", $8, $10, $12;
	else
	    print ls, ev, $8, $10, $12, " ", prob, mlb, mcmatch;
    };
    ls=$4;
    ev=$6;
    prob=$8;
    mlb=$10;
}

