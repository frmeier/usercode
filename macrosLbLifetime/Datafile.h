#ifndef DATAFILE_H_GUARD
#define DATAFILE_H_GUARD

#include <string>
#include <vector>
#include "TChain.h"
#include "HtmlReport.h"
#include "utils.h"

struct Datafile
{
    public:
	Datafile(std::string setfilename, std::string settitle, double setlumi, double setevents, int setrunMin, int setrunMax) :
	    filename(setfilename), title(settitle), lumi(setlumi), events(setevents), runMin(setrunMin), runMax(setrunMax) {};

	std::string filename;
	std::string title;
	double lumi; // in /pb
	double events;
	int runMin;
	int runMax;
};

struct Datafiles
{
    public:
	void add(Datafile d) { datafilevec.push_back(d); };
	double getLumi()
	{
	    double ret(0);
	    for (std::vector<Datafile>::const_iterator it = datafilevec.begin(); it!= datafilevec.end(); it++)
	    {
		ret += it->lumi;
	    }
	    return ret;
	};
	double getLumiPb() { return getLumi(); };
	double getLumiPbRounded() { return int(getLumiPb()+.5); };
	double getEvents()
	{
	    double ret(0);
	    for (std::vector<Datafile>::const_iterator it = datafilevec.begin(); it!= datafilevec.end(); it++)
	    {
		ret += it->events;
	    }
	    return ret;
	};
	TChain* getTChain(std::string treename)
	{
	    TChain* ret = new TChain(treename.c_str());
	    for (std::vector<Datafile>::const_iterator it = datafilevec.begin(); it!= datafilevec.end(); it++)
	    {
		ret->Add(it->filename.c_str());
	    }
	    return ret;
	}
	void mkHtmlReport(HtmlReport *ht, std::string httitle)
	{
	    ht->addH(httitle,'3');
	    ht->addTableCell("File");
	    ht->addTableCell("Lumi pb-1");
	    ht->addTableCell("Events");
	    ht->endTableRow();
	    for (std::vector<Datafile>::const_iterator it = datafilevec.begin(); it != datafilevec.end(); it++)
	    {
		ht->addTableRow(it->filename, roundToString(it->lumi*0.000001,1), toString(it->events));
	    }
	    ht->addTableRow("Total:", toString(getLumiPbRounded()), toString(getEvents()));
	};

	std::string title;
	std::vector<Datafile> datafilevec;
};

#endif
