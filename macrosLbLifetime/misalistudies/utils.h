#ifndef UTILS_H_GUARD
#define UTILS_H_GUARD

#include "TLatex.h"
#include <sstream>
#include <iostream>
#include <string>
#include <iostream>
#include <iomanip>
#include <vector>

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

std::string makePlotsString(const std::string &toPlot, const std::string &hname, const int &nBins, const double &lo, const double &hi)
{
    return toPlot + ">>" + hname + "(" + toString(nBins) + "," + toString(lo) + "," + toString(hi) + ")";
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

#endif
