#include <string>
#include <vector>
#include <utility>
#include <iostream>

#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1F.h"
#include "TEventList.h"
#include "TPad.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLorentzVector.h"

#include "utils.h"
#include "Cuts.C"
#include "do1dPlot.C"
#include "do1dPlotGaus.C"
#include "setTDRStyle_modified.C"

void doResoPlots(string filename, string title)
{
    setTDRStyle();
    TFile *f = TFile::Open(filename.c_str()); // Data
    TTree *events = (TTree*)f->Get("events");

    Cuts cut;
    //cut.selectCut("lb12", "acc06Lb", "muSoft", "HLT_jpsiBarrel", "HLT_matched");
    cut.selectCut("B006", "acc06B0", "muSoft", "HLT_jpsiBarrel", "HLT_matched");
    //const string curCut = cut.getCut() + "&&mbc>5.24&&mbc<5.32";
    //const string curCut = "isMCmatch==1&&mbc>5.24&&mbc<5.32";
    const string curCut = "isMCmatch==1&&rqha1<0";
    //const string curCut = cut.getCut();

    //const double ctWindow(10e-12), ctWindowZoom(1e-12);
    const double ctWindow(10), ctWindowZoom(1);
    //do1dPlot(events, "MCtruthCt3dLb", "ct3dlb-ctlbtruth", "", 40, -ctWindow, ctWindow, "Lifetime of #Lambda_{b} 3d, all candidates found", "#tau_{reco,3d}(#Lambda_{b})-#tau_{truth}(#Lambda_{b})", "s");
    //do1dPlot(events, "MCtruthCt3dLb", "ct3dbc-ctbctruth", "", 40, -ctWindowZoom, ctWindowZoom, "Lifetime of B^{0} 3d, all candidates found", "#tau_{reco,3d}(B^{0})-#tau_{truth}(B^{0})", "s");
    //do1dPlotGaus(events, "MCtruthCt3dLb", "ct3dbc-ctbctruth", "", 40, -ctWindowZoom, ctWindowZoom, "Lifetime of B^{0} 3d, all candidates found", "#tau_{reco,3d}(B^{0})-#tau_{truth}(B^{0})", "s");

    do1dPlotGaus(events, "MCtruthCt3dLb", "1e12*(ct3dbc-ctbctruth)", "", 40, -ctWindowZoom, ctWindowZoom, title, "#tau_{reco,3d}(B^{0})-#tau_{truth}(B^{0})", "ps");
    //do1dPlotGaus(events, "MCtruthCt3dLb", "1e12*(ct3dbc-ctbctruth)", "", 40, -ctWindowZoom, ctWindowZoom, title, "#tau_{reco,3d}(#Lambda_{b})-#tau_{truth}(#Lambda_{b})", "ps");
    //do1dPlotGaus(events, "MCtruthCt3dLb", "1e12*(ct3dbc-ctbctruth)", "rqha1<0", 40, -ctWindowZoom, ctWindowZoom, title, "#tau_{reco,3d}(#Lambda_{b})-#tau_{truth}(#Lambda_{b})", "ps");
}

void doSomeResoPlots()
{
    const string path = "../data/";
    std::vector< std::pair<string,string> > filelist;

    filelist.push_back(std::make_pair("radialEpsilon_minus", path+"run491.root"));
    filelist.push_back(std::make_pair("radialEpsilon_plus", path+"run492.root"));
    filelist.push_back(std::make_pair("telescopeEpsilon_minus", path+"run497.root"));
    filelist.push_back(std::make_pair("telescopeEpsilon_plus", path+"run498.root"));
    filelist.push_back(std::make_pair("layerRotEpsilon_minus", path+"run489.root"));
    filelist.push_back(std::make_pair("layerRotEpsilon_plus", path+"run490.root"));

    filelist.push_back(std::make_pair("bowingEpsilon_minus", path+"run485.root"));
    filelist.push_back(std::make_pair("bowingEpsilon_plus", path+"run486.root"));
    filelist.push_back(std::make_pair("zExpEpsilon_minus", path+"run501.root"));
    filelist.push_back(std::make_pair("zExpEpsilon_plus", path+"run502.root"));
    filelist.push_back(std::make_pair("twistEpsilon_minus", path+"run499.root"));
    filelist.push_back(std::make_pair("twistEpsilon_plus", path+"run500.root"));

    filelist.push_back(std::make_pair("ellipticalEpsilon_minus", path+"run487.root"));
    filelist.push_back(std::make_pair("ellipticalEpsilon_plus", path+"run488.root"));
    filelist.push_back(std::make_pair("skewEpsilon_minus2", path+"run527.root"));
    filelist.push_back(std::make_pair("skewEpsilon_plus2", path+"run528.root"));
    filelist.push_back(std::make_pair("saggitaEpsilon_minus", path+"run493.root"));
    filelist.push_back(std::make_pair("saggitaEpsilon_plus", path+"run494.root"));

    /*
    filelist.push_back(std::make_pair("skewEpsilon_minus", path+"run495.root"));
    filelist.push_back(std::make_pair("skewEpsilon_plus", path+"run496.root"));
    filelist.push_back(std::make_pair("radialEpsilon_layer1minus", path+"run530.root"));
    filelist.push_back(std::make_pair("radialEpsilon_layer1plus", path+"run531.root"));
    filelist.push_back(std::make_pair("surfdeform", path+"run529.root"));
    */

    /*
    filelist.push_back(std::make_pair("radialEpsilon_minus", path+"run509.root"));
    filelist.push_back(std::make_pair("radialEpsilon_plus", path+"run510.root"));
    filelist.push_back(std::make_pair("telescopeEpsilon_minus", path+"run515.root"));
    filelist.push_back(std::make_pair("telescopeEpsilon_plus", path+"run516.root"));
    filelist.push_back(std::make_pair("layerRotEpsilon_minus", path+"run507.root"));
    filelist.push_back(std::make_pair("layerRotEpsilon_plus", path+"run508.root"));

    filelist.push_back(std::make_pair("bowingEpsilon_minus", path+"run503.root"));
    filelist.push_back(std::make_pair("bowingEpsilon_plus", path+"run504.root"));
    filelist.push_back(std::make_pair("zExpEpsilon_minus", path+"run519.root"));
    filelist.push_back(std::make_pair("zExpEpsilon_plus", path+"run520.root"));
    filelist.push_back(std::make_pair("twistEpsilon_minus", path+"run517.root"));
    filelist.push_back(std::make_pair("twistEpsilon_plus", path+"run518.root"));

    filelist.push_back(std::make_pair("ellipticalEpsilon_minus", path+"run505.root"));
    filelist.push_back(std::make_pair("ellipticalEpsilon_plus", path+"run506.root"));
    filelist.push_back(std::make_pair("skewEpsilon_minus2", path+"run522.root"));
    filelist.push_back(std::make_pair("skewEpsilon_plus2", path+"run523.root"));
    filelist.push_back(std::make_pair("saggitaEpsilon_minus", path+"run511.root"));
    filelist.push_back(std::make_pair("saggitaEpsilon_plus", path+"run512.root"));
    */

    /*
    filelist.push_back(std::make_pair("skewEpsilon_minus", path+"run513.root"));
    filelist.push_back(std::make_pair("skewEpsilon_plus", path+"run514.root"));
    filelist.push_back(std::make_pair("radialEpsilon_layer1minus", path+"run525.root"));
    filelist.push_back(std::make_pair("radialEpsilon_layer1plus", path+"run526.root"));
    filelist.push_back(std::make_pair("surfdeform", path+"run524.root"));
    */

    TCanvas *c = new TCanvas("c", "c", 800, 800);
    c->Divide(6,3);

    for (unsigned int i = 0; i!=filelist.size(); i++)
    {
	cout << "---- " << filelist[i].first << endl;
	c->cd(i+1);
	doResoPlots(filelist[i].second, filelist[i].first);
    }
}

/*
bowingEpsilon_minus
bowingEpsilon_plus
ellipticalEpsilon_minus
ellipticalEpsilon_plus
layerRotEpsilon_minus
layerRotEpsilon_plus
radialEpsilon_minus
radialEpsilon_plus
saggitaEpsilon_minus
saggitaEpsilon_plus
skewEpsilon_minus
skewEpsilon_plus
telescopeEpsilon_minus
telescopeEpsilon_plus
twistEpsilon_minus
twistEpsilon_plus
zExpEpsilon_minus
zExpEpsilon_plus
*/

