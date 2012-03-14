#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "THStack.h"
#include "TLatex.h" 
#include "setTDRStyle_modified.C"
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <iomanip>

using std::cout;
using std::endl;

TCanvas *c1;

struct strctHistoStacked{
    strctHistoStacked(std::string pFilename, std::string pName, std::string pTitle, int pColor, std::string pCut)
    {
	valid = false;
	file = TFile::Open(pFilename.c_str());
	title = pTitle;
	name = pName;
	color = pColor;
	cut = pCut;
	if (file != 0)
	{
	    valid = true;
	    tree = (TTree*) file->Get("events");
	}
    };
    std::string title;
    std::string name;
    std::string cut;
    int color;
    TFile* file;
    TTree* tree;
    TH1F* histo;
    TPaveStats* statbox;
    bool valid;
};

typedef std::vector<strctHistoStacked> strctHistoStacked_type;

template <typename T>
std::string toString(T i)
{
    ostringstream oss;
    oss << i;
    return oss.str();
}

std::string formatString(double val, std::streamsize prec)
{
    ostringstream oss;
    std::streamsize oldprec = oss.precision();
    oss << std::setprecision(prec) << val << std::setprecision(oldprec);
    return oss.str();
}

void setStageStacked(std::string strPlot, std::string strCut, int nBins, double mymin, double mymax, bool stacked = true, bool doNorm = false) {
    setTDRStyle();

    gStyle->SetOptStat(112211);

    // some general settings for the plot
    //std::string mainTitle = "Run range study - MuOnia L=19.6/pb";
    std::string mainTitle = "L1 trigger in MuOnia";
    std::string xAxisUnit = "GeV/c^{2}";
    //std::string xAxisTitle = strPlot + " / " + xAxisUnit;
    std::string xAxisTitle = "m(#mu#mup#pi) / " + xAxisUnit;
    double binsize = (mymax-mymin)/nBins;
    std::string strBinsize = formatString(binsize,3);

    // Set up vector with histos
    strctHistoStacked_type vecHistos;
    std::string strPath = "/scratch/frmeier/";
    //std::string strPath = "dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/frmeier/readerRuns/mergedTrees/";
    //std::string strPath = "dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/frmeier/readerRuns/run059/";

    /*
    vecHistos.push_back(strctHistoStacked::strctHistoStacked((strPath+"run026.root").c_str(),"0","443 3122",1,"nRef1G==0"));
    vecHistos.push_back(strctHistoStacked::strctHistoStacked((strPath+"run026.root").c_str(),"1","443 3122 211 211",2,"nRef1G==1"));
    vecHistos.push_back(strctHistoStacked::strctHistoStacked((strPath+"run026.root").c_str(),"2","443 3122 211 211 22",3,"nRef1G==2"));
    vecHistos.push_back(strctHistoStacked::strctHistoStacked((strPath+"run026.root").c_str(),"3","443 3122 211 211 22 22",4,"nRef1G==3"));
    vecHistos.push_back(strctHistoStacked::strctHistoStacked((strPath+"run026.root").c_str(),"gt3","others",6,"nRef1G>3"));
    */

    //vecHistos.push_back(strctHistoStacked::strctHistoStacked((strPath+"run047.root").c_str(),"hsig047_e","#Lambda_{b}#rightarrowJ/#psi(#mu#mu)#Lambda(p#piX)",38,"nRef2G>=921&&nRef2G<=923"));
    //vecHistos.push_back(strctHistoStacked::strctHistoStacked((strPath+"run047.root").c_str(),"hsig047_f","#Lambda_{b}#rightarrowJ/#psi(#mu#mu)#Lambda(others)",20,"((nRef2G>=216&&nRef2G<511)||(nRef2G>923&&nRef2G<=1254))"));

    /*
    vecHistos.push_back(strctHistoStacked::strctHistoStacked((strPath+"run047.root").c_str(),"hsig047_a","#Lambda_{b}#rightarrow other",45,"(nRef2G<216||nRef2G>4811)"));
    vecHistos.push_back(strctHistoStacked::strctHistoStacked((strPath+"run047.root").c_str(),"hsig047_b","#Lambda_{b}#rightarrowJ/#psi(#mu#mu)X(other)",41,"nRef2G>1254&&nRef2G<=4811"));
    vecHistos.push_back(strctHistoStacked::strctHistoStacked((strPath+"run047.root").c_str(),"hsig047_c","#Lambda_{b}#rightarrowJ/#psi(#mu#mu)#Lambda(p#pi)X",38,"nRef2G>511&&nRef2G<=920"));
    vecHistos.push_back(strctHistoStacked::strctHistoStacked((strPath+"run047.root").c_str(),"hsig047","#Lambda_{b}#rightarrowJ/#psi(#mu#mu)#Lambda(p#pi)",3,"nRef2G==511"));
    */

    /*
    vecHistos.push_back(strctHistoStacked::strctHistoStacked((strPath+"run052.root").c_str(),"hsig052_a","isSig==0&&isMCmatch==0",1,"isSig==0&&isMCmatch==0"));
    vecHistos.push_back(strctHistoStacked::strctHistoStacked((strPath+"run052.root").c_str(),"hsig052_c","isSig==0&&isMCmatch==1",4,"isSig==0&&isMCmatch==1"));
    vecHistos.push_back(strctHistoStacked::strctHistoStacked((strPath+"run052.root").c_str(),"hsig052_d","isSig==1&&isMCmatch==1",3,"isSig==1&&isMCmatch==1"));
    vecHistos.push_back(strctHistoStacked::strctHistoStacked((strPath+"run052.root").c_str(),"hsig052_b","isSig==1&&isMCmatch==0",2,"isSig==1&&isMCmatch==0"));
    */

    //vecHistos.push_back(strctHistoStacked::strctHistoStacked((strPath+"run047.root").c_str(),"hdat039_hi","second half",3,"run>=147283"));

    //vecHistos.push_back(strctHistoStacked::strctHistoStacked((strPath+"odd.chain.root").c_str(),"hdat058odd","odd",2,""));
    //vecHistos.push_back(strctHistoStacked::strctHistoStacked((strPath+"even.chain.root").c_str(),"hdat058even","even",3,""));

    // Triggerstudie
    vecHistos.push_back(strctHistoStacked::strctHistoStacked((strPath+"run090.root").c_str(),"hdat090","alle",1,""));
    vecHistos.push_back(strctHistoStacked::strctHistoStacked((strPath+"run090.root").c_str(),"hdat090L1SMu0DMu0","SMu0||DMu0",7,"((L1TMu0==1)||(L1TDMu0==1))"));
    vecHistos.push_back(strctHistoStacked::strctHistoStacked((strPath+"run090.root").c_str(),"hdat090L1TMu0","L1_SingleMu0",6,"L1TMu0==1"));
    vecHistos.push_back(strctHistoStacked::strctHistoStacked((strPath+"run090.root").c_str(),"hdat090L1TDMu0","L1_DoubleMuOpen",2,"L1TDMu0==1"));
    vecHistos.push_back(strctHistoStacked::strctHistoStacked((strPath+"run090.root").c_str(),"hdat090L1TDMu3","L1_DoubleMu3",3,"L1TDMu3==1"));
    vecHistos.push_back(strctHistoStacked::strctHistoStacked((strPath+"run090.root").c_str(),"hdat090L1TMu10","L1_SingleMu10",7,"L1TMu10==1"));
    //vecHistos.push_back(strctHistoStacked::strctHistoStacked((strPath+"run090.root").c_str(),"hdat090","",3,"==1"));
    //vecHistos.push_back(strctHistoStacked::strctHistoStacked((strPath+"run090.root").c_str(),"hdat090","",3,"==1"));
    //vecHistos.push_back(strctHistoStacked::strctHistoStacked((strPath+"run090.root").c_str(),"hdat090","",3,"==1"));
    //vecHistos.push_back(strctHistoStacked::strctHistoStacked((strPath+"run090.root").c_str(),"hdat090","",3,"==1"));
    //vecHistos.push_back(strctHistoStacked::strctHistoStacked((strPath+"run090.root").c_str(),"hdat090","",3,"==1"));


    // check if the files were opened correctly
    {
	strctHistoStacked_type::iterator it = vecHistos.begin();
	while (it != vecHistos.end())
	{
	    if (it->valid)
	    {
		it++;
	    }
	    else
	    {
		cout << "The file for the histo " << it->name
		    << " " << it->title
		    << " could not be opened - histo dropped from list."
		    << endl;
		it = vecHistos.erase(it);
	    }
	}
    }

    // check no of histos and complain if problems may arise
    int nHistos = vecHistos.size();
    if (nHistos > 5) {
	cout << "Warning: You request to superimpose more than 5 histos." << endl
	   << "I'm running out of space for the statboxes..." << endl;
    }
    if (nHistos == 0)
    {
	cout << "No histograms in list. Ending." << endl;
	return;
    }

    // Canvas
    c1 = new TCanvas;
    c1->SetWindowSize(1200,600);
    TPad* pad1= (TPad*)c1->cd();
    pad1->cd();
    pad1->SetTicks(0,0);
    pad1->SetRightMargin(0.20);
    pad1->SetTopMargin(0.10);

    // palette for colors
    //const int colorsarr[] = {1,2,3,4,6,7,8,9};
    const int colorsarr[] = {1,4,9,38,32,30,8};
    int icolor(0);
    const bool useColorsarr(true);
    const int hatcharr[] = {3515,3535,3545};

    // plotting
    std::string strHist("("+toString(nBins)+","+toString(mymin)+","+toString(mymax)+")");
    cout << strHist << endl;
    THStack *hstck = new THStack("hs2",mainTitle.c_str());
    for (strctHistoStacked_type::iterator it = vecHistos.begin(); it!= vecHistos.end(); it++, icolor++)
    {
	std::string curCut;
	if ((*it).cut.size() != 0)
	    curCut = strCut + "&&" + (*it).cut;
	else
	    curCut = strCut;
	cout << curCut << endl;
	// draw, but the first one needs no "sames"
	(*it).tree->Draw((strPlot+">>"+(*it).name+strHist).c_str(), curCut.c_str());
	(*it).histo = (TH1F*)gDirectory->GetList()->FindObject((*it).name.c_str());
	(*it).histo->SetLineColor(1);
	//(*it).histo->SetFillStyle(hatcharr[icolor % (sizeof(hatcharr)/sizeof(int))]);
	if(useColorsarr)
	    (*it).histo->SetFillColor(colorsarr[icolor % (sizeof(colorsarr)/sizeof(int))]);
	else
	    (*it).histo->SetFillColor((*it).color);
	//(*it).histo->GetXaxis()->SetTitle(xAxisTitle.c_str());

	//(*it).histo->GetYaxis()->SetTitle(("events per " + strBinsize + " " + xAxisUnit).c_str());
	// normalize if required
	if(doNorm) (*it).histo->Scale(1./(*it).histo->GetSumOfWeights());
	hstck->Add((*it).histo,"sames");
    }
    
    // draw the histos unstacked to get the correct statboxes
    hstck->Draw("HIST nostack");

    // axis title
    hstck->GetXaxis()->SetTitle(xAxisTitle.c_str());
    hstck->GetXaxis()->SetTicks("+-");
    hstck->GetXaxis()->SetTitleOffset(1.0);
    hstck->GetYaxis()->SetTitle(("events per " + strBinsize + " " + xAxisUnit).c_str());
    hstck->GetYaxis()->SetTicks("+-");
    hstck->GetYaxis()->SetTickLength(0.01);
    hstck->GetYaxis()->SetTitleOffset(-0.8);

    // Statboxen
    pad1->Update();
    double deltay=0, top_corner = 0.9, bottom_corner = .13;
    double dHeightStatbox = .6/((double)nHistos);
    TPaveStats *st1;
    if (stacked)
    {
	icolor=vecHistos.size()-1;
	for (strctHistoStacked_type::reverse_iterator it = vecHistos.rbegin(); it!= vecHistos.rend(); it++, icolor--)
	{	// we do this in reverse order as this will give a more natural sequence for a stacked histo
	    st1 = (TPaveStats*)(*it).histo->GetListOfFunctions()->FindObject("stats");
	    st1->SetOptStat(112211);
	    st1->SetY1NDC(top_corner-deltay-dHeightStatbox);
	    st1->SetY2NDC(top_corner-deltay);
	    st1->SetX1NDC(0.80); st1->SetX2NDC(.995);
	    if(useColorsarr)
		st1->SetTextColor(colorsarr[icolor % (sizeof(colorsarr)/sizeof(int))]);
	    else
		st1->SetTextColor((*it).color);
	    //st1->SetStatFormat("2.3e");
	    //st1->SetFitFormat("2.3e");
	    deltay+=dHeightStatbox;
	    (*it).statbox = (TPaveStats*)st1->Clone();
	}
    }
    else
    {
	icolor=0;
	for (strctHistoStacked_type::iterator it = vecHistos.begin(); it!= vecHistos.end(); it++, icolor++)
	{	// we do this in reverse order as this will give a more natural sequence for a stacked histo
	    st1 = (TPaveStats*)(*it).histo->GetListOfFunctions()->FindObject("stats");
	    st1->SetOptStat(112211);
	    st1->SetY1NDC(top_corner-deltay-dHeightStatbox);
	    st1->SetY2NDC(top_corner-deltay);
	    st1->SetX1NDC(0.80); st1->SetX2NDC(.995);
	    if(useColorsarr)
		st1->SetTextColor(colorsarr[icolor % (sizeof(colorsarr)/sizeof(int))]);
	    else
		st1->SetTextColor((*it).color);
	    //st1->SetStatFormat("2.3e");
	    //st1->SetFitFormat("2.3e");
	    deltay+=dHeightStatbox;
	    (*it).statbox = (TPaveStats*)st1->Clone();
	}
    }

    // Legende einfügen
    TLegend *legend = new TLegend(0.80,top_corner-deltay,0.995,bottom_corner);
    legend->SetFillStyle(1000);
    legend->SetBorderSize(1.);
    legend->SetTextSize(0.02);
    legend->SetFillColor(0);
    //legend->SetBorderSize(0.);
    if (stacked)
    {
	for (strctHistoStacked_type::reverse_iterator it = vecHistos.rbegin(); it!= vecHistos.rend(); it++)
	{
	    legend->AddEntry((*it).histo,(*it).title.c_str(),"fl");
	}
    }
    else
    {
	for (strctHistoStacked_type::iterator it = vecHistos.begin(); it!= vecHistos.end(); it++)
	{
	    legend->AddEntry((*it).histo,(*it).title.c_str(),"fl");
	}
    }


    // draw
    pad1->Update();
    pad1->Modified();
    pad1->Draw("");
    if (stacked)
	hstck->Draw("HIST");
    else
	hstck->Draw("HIST nostack");
    for (strctHistoStacked_type::iterator it = vecHistos.begin(); it!= vecHistos.end(); it++)
    {
	(*it).statbox->Draw();
    }
    legend->Draw();

    gPad->RedrawAxis();
    //pad1->Update();

    // Text for stacked or not stacked
    TLatex txt;
    txt.SetTextSize(0.02);
    txt.SetTextAlign(13);
    txt.SetNDC(true);
    txt.DrawLatex(0.807,bottom_corner-.01,stacked ? "stacked":"overlayed, not stacked");

    c1->SaveAs("c1.pdf");
}

