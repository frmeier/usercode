BEGIN {
    v=0.2;
    vinc=.2;
    bg = 0;
    bge = 0;
}

/Lifetime/ {
    printf "%s %s %s ", v, $2, $4;
    getline;
    printf "%s %s %s  bg: %s %s\n", $1, $2, $3, bg, bge;
    v+=vinc;
}

/12  tau_bk/ {
    bg = $3*1e12;
    bge= $4*1e12;
}

