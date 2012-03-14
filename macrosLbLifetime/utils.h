#ifndef UTILS_H_GUARD
#define UTILS_H_GUARD

#include "TLatex.h"
#include <sstream>
#include <iostream>
#include <string>
#include <iostream>
#include <iomanip>

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

#endif
