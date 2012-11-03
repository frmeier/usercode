#ifndef UTILS_H_GUARD
#define UTILS_H_GUARD

#include "TROOT.h"
#include "TLatex.h"
#include "TTree.h"
#include "TEventList.h"
#include <sstream>
#include <iostream>
#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include <list>

template <typename T>
std::string toString(T i)
{
    std::ostringstream oss;
    oss << i;
    return oss.str();
};

double square(double v) { return v*v; };

std::string roundToString(double v, std::streamsize precision)
{
    std::ostringstream oss;
    oss << std::setprecision(precision) << std::fixed << v;
    return oss.str();
}

TLatex* writeTLatex(std::string text, double x, double y, double size = 0.04)
{
    TLatex* txt = new TLatex(x,y,text.c_str());
    txt->SetTextSize(size);
    txt->SetNDC(true);
    return txt;
}

std::string muIdStr(const int &id)
{
    const std::string sid = toString(id);
    return "(rid1m&"+sid+")=="+sid+"&&(rid2m&"+sid+")=="+sid;
}

std::string valueWithUnit(const std::string value, const std::string unit)
{
    if (unit.size() != 0)
	return value + " [" + unit + "]";
    else
	return value;
}

std::string entriesPerBin(const double binsize, std::string unit)
{
    return "entries per " + valueWithUnit(toString(binsize), unit);
}

std::string entriesPerBin(const int nBins, const double min, const double max, std::string unit)
{
    const double binsize = (max-min)/nBins;
    return entriesPerBin(binsize, unit);
}

std::string concatCutString(const std::string s1, const std::string s2, const std::string s3 = "", const std::string s4 = "", const std::string s5 = "")
{
    const std::string ret = s1 + "&&" + s2;
    if (s3.size() != 0)
	return concatCutString(ret, s3, s4, s5);
    else
	return ret;
}

std::string makePlotsString(const int &nBins, const double &lo, const double &hi)
{
    return "(" + toString(nBins) + "," + toString(lo) + "," + toString(hi) + ")";
}

std::string makePlotsString(const std::string &toPlot, const std::string &hname, const int &nBins, const double &lo, const double &hi)
{
    return toPlot + ">>" + hname + "(" + toString(nBins) + "," + toString(lo) + "," + toString(hi) + ")";
}

std::string makePlotsString(const std::string &toPlot, const std::string &hname)
{
    return toPlot + ">>" + hname;
}

std::vector<double> variableBinSizeVec(double v0, double v1, double v2 = -9999, double v3 = -9999, double v4 = -9999, double v5 = -9999, double v6 = -9999,
	double v7 = -9999, double v8 = -9999, double v9 = -9999, double v10 = -9999, double v11 = -9999, double v12 = -9999, double v13 = -9999)
{
    std::vector<double> ret;
    ret.push_back(v0);
    ret.push_back(v1);
    if (v2 != -9999) ret.push_back(v2); else return ret;
    if (v3 != -9999) ret.push_back(v3); else return ret;
    if (v4 != -9999) ret.push_back(v4); else return ret;
    if (v5 != -9999) ret.push_back(v5); else return ret;
    if (v6 != -9999) ret.push_back(v6); else return ret;
    if (v7 != -9999) ret.push_back(v7); else return ret;
    if (v8 != -9999) ret.push_back(v8); else return ret;
    if (v9 != -9999) ret.push_back(v9); else return ret;
    if (v10 != -9999) ret.push_back(v10); else return ret;
    if (v11 != -9999) ret.push_back(v11); else return ret;
    if (v12 != -9999) ret.push_back(v12); else return ret;
    if (v13 != -9999) ret.push_back(v13); else return ret;
    return ret;
}

std::vector<double> makeDynamicBins(TTree *t, string branchname, string cut, int nBins, double lo, double hi, const double factor = 1)
{
    typedef std::list<double> listtype;
    listtype valuelist;

    t->Draw(">>lst", cut.c_str());
    TEventList *lst = (TEventList*)gDirectory->Get("lst");

    double value;
    t->SetBranchAddress(branchname.c_str(), &value);
    const int N = lst->GetN();

    for (int i=0; i!=N; i++)
    {
	t->GetEntry(lst->GetEntry(i));
	value*=factor;
	if (value>=lo && value <=hi)
	{
	    valuelist.push_back(value);
	}
    }
    valuelist.sort();
    const int entriesPerBin = valuelist.size() / nBins;

    std::vector<double> binvector;
    binvector.push_back(lo);

    int i = 0;
    for (listtype::const_iterator it = valuelist.begin(); it!=valuelist.end(); it++, i++)
    {
	if (i==0) continue;
	if (i % entriesPerBin == 0) binvector.push_back(*it);
    }
    binvector.pop_back(); // remove last one (ugly, but easier than awful logic in loop)
    binvector.push_back(hi); // and replace it with given boundary
    //for (std::vector<double>::const_iterator it = binvector.begin(); it!=binvector.end(); it++)
    //	cout << *it << " ";
    //cout << endl;
    delete lst;
    t->ResetBranchAddresses();

    return binvector;
}

#endif
