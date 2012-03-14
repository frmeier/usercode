BEGIN {
    n=0;
}

/Using additional cut/ { 
    for (i=0; i!=n; i++) { print  a[i]; };
    print "-------------------------------------------------------------------";
    print;
} 

/curEntries:/ { print }

/EXT PARAMETER/ { n=0;}

/EXT PARAMETER/,/ERR DEF/ { 
    a[n]=$0;
    #print n ": " $0;
    n++;
}
    
END {
    for (i=0; i!=n; i++) { print  a[i]; };
}

