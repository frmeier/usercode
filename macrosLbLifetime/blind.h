#ifndef GUARD_BLIND_H
#define GUARD_BLIND_H

#include "TRandom3.h"

double blindFactor()
{
    TRandom3 ran;
    ran.SetSeed(17+18*3-5+17);
    return ran.Uniform(0.8,1.2);
}

#endif
