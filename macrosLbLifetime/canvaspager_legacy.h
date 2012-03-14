#ifndef CANVASPAGER_LEGACY_H_GUARD
#define CANVASPAGER_LEGACY_H_GUARD
#include "utils.h"

void canvaspager(TCanvas *c, std::string filename, const int &ncd, int &cdcounter, int &pagecounter)
{
    cdcounter++;
    if (cdcounter > ncd)
    {
	c->SaveAs((filename + toString(pagecounter) + ".pdf").c_str());
	c->SaveAs((filename + toString(pagecounter) + ".png").c_str());
	pagecounter++;
	cdcounter=1;
	c->Clear("D");
    }
}

#endif

