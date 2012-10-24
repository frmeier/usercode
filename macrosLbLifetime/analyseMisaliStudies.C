#include <string>
#include <iostream>
#include <vector>
#include <utility>

#include "Fit_2D.C"

void analyseMisaliStudies()
{
    typedef vector< pair<string,string> > list_t;
    const string path = "../data/";
    list_t liste;
    //liste.push_back(make_pair("569","bowingEpsilon_minus"));
    //liste.push_back(make_pair("570","bowingEpsilon_plus"));
    //liste.push_back(make_pair("571","ellipticalEpsilon_minus"));
    liste.push_back(make_pair("572","ellipticalEpsilon_plus"));
    liste.push_back(make_pair("573","layerRotEpsilon_minus"));
    liste.push_back(make_pair("574","layerRotEpsilon_plus"));
    liste.push_back(make_pair("575","nomisali"));
    //liste.push_back(make_pair("576","radialEpsilon_layer1minus2"));
    liste.push_back(make_pair("577","radialEpsilon_layer1plus2"));
    liste.push_back(make_pair("578","radialEpsilon_minus"));
    liste.push_back(make_pair("579","radialEpsilon_plus"));
    liste.push_back(make_pair("580","saggitaEpsilon_minus"));
    liste.push_back(make_pair("581","saggitaEpsilon_plus"));
    //liste.push_back(make_pair("582","skewEpsilon_minus2"));
    //liste.push_back(make_pair("583","skewEpsilon_plus2"));
    liste.push_back(make_pair("584","telescopeEpsilon_minus"));
    liste.push_back(make_pair("585","telescopeEpsilon_plus"));
    liste.push_back(make_pair("586","twistEpsilon_minus"));
    liste.push_back(make_pair("587","twistEpsilon_plus"));
    liste.push_back(make_pair("588","zExpEpsilon_minus"));
    liste.push_back(make_pair("589","zExpEpsilon_plus"));
    
    for (list_t::const_iterator it=liste.begin(); it!=liste.end(); it++)
    {
	TString filename(path+"vrt_r"+it->first+"_"+it->second+".root");
	Fit_2D(filename);
    }
}

