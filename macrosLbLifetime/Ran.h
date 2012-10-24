#include "TRandom3.h"

TRandom3 rn;
double getrn() { return rn.Uniform(); }

double getProb(double r)
{
    if (r<3) return .3+.117*r;
    if (r<6) return .65;
    if (r<7) return .65+(r-6)*.06;
    if (r<8) return .71-(r-7)*.13;
    if (r<12) return .58;
    if (r<24) return .58+(r-12)*.035;
    if (r<35) return 1-(r-24)*.038;
    if (r<50) return .58;
    if (r<80) return .58-(r-50)*.02;
    return 0.0;
}

double getProbZ(double z)
{
    z = TMath::Abs(z);
    if (z<5) return 1.0;
    if (z<40) return 1-(z-5)*.01;
    if (z<80) return .667;
    if (z<90) return .667+(z-80)*.017;
    if (z<200) return .84-(z-90)*.008;
    return 0.0;
}

double getProbR(double r)
{
    if (r<.1)
	return .02*.08*r;
    else
	return 1;
}

/*

class Ran
{
    public:
	double get() { return ran_.Uniform(); };

    private:
	TRandom3 ran_;
};

#if defined(__MAKECINT__)
#pragma link C++ class Ran;
#endif
*/

