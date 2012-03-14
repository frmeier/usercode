#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "setTDRStyle_modified.C"
#include <string>
#include <sstream>
#include <iostream>
#include <vector>

using std::cout;
using std::endl;

TCanvas *c1;

struct strctHisto{
    strctHisto(std::string pFilename, std::string pName, std::string pTitle, int pColor)
    {
	valid = false;
	file = TFile::Open(pFilename.c_str());
	title = pTitle;
	name = pName;
	color = pColor;
	if (file != 0)
	{
	    valid = true;
	    tree = (TTree*) file->Get("events");
	}
    };
    std::string title;
    std::string name;
    int color;
    TFile* file;
    TTree* tree;
    TH1F* histo;
    bool valid;
};

typedef std::vector<strctHisto> strctHisto_type;

template <typename T>
std::string toString(T i)
{
    ostringstream oss;
    oss << i;
    return oss.str();
}

void setStage(std::string strPlot, std::string strCut, int nBins, double mymin, double mymax, bool doNorm = false) {
    setTDRStyle();

    gStyle->SetOptStat(112211);

    // Colors
    std::vector<int> colors(3,1);
    colors[0]=2; colors[1]=4;

    // Set up vector with histos
    strctHisto_type vecHistos;
    std::string strPath = "/scratch/frmeier/";
    //std::string strPath = "/media/DR-1/";
    //std::string strPath = "/home/frank/psi/lambda/data/";
    //std::string strPath = "dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/frmeier/readerRuns/mergedTrees/";
    //std::string strPath = "/shome/meier_f1/CMSSW/CMSSW_3_8_4_patch2/src/HeavyFlavorAnalysis/Bs2MuMu/macros/";

    /*
    vecHistos.push_back(strctHisto::strctHisto((strPath+"run009.root").c_str(),"hsig6","Signal MC Cand 6",1));
    vecHistos.push_back(strctHisto::strctHisto((strPath+"run012.root").c_str(),"hsig7","Signal MC Cand 7",2));
    vecHistos.push_back(strctHisto::strctHisto((strPath+"run015.root").c_str(),"hsig8","Signal MC Cand 8",3));
    vecHistos.push_back(strctHisto::strctHisto((strPath+"run018.root").c_str(),"hsig9","Signal MC Cand 9",4));
    */
    //vecHistos.push_back(strctHisto::strctHisto((strPath+"run009part.root").c_str(),"hsig006","Signal MC C6",1));

    //vecHistos.push_back(strctHisto::strctHisto("/scratch/frmeier/run011part.root","hbgr6","Background MC",2));
    //vecHistos.push_back(strctHisto::strctHisto("/scratch/frmeier/run010.root","hdat6","Data",1));
    //vecHistos.push_back(strctHisto::strctHisto((strPath+"run019.root").c_str(),"hdat019","Data old",4));
    //vecHistos.push_back(strctHisto::strctHisto((strPath+"run031.root").c_str(),"hdat031","Data new",2));

    //vecHistos.push_back(strctHisto::strctHisto((strPath+"run033.root").c_str(),"hdat033","MuOnia, L=12.2/pb",1));
    //vecHistos.push_back(strctHisto::strctHisto((strPath+"run034.root").c_str(),"hdat034","MuOnia, L=14.3/pb",4));
    //vecHistos.push_back(strctHisto::strctHisto((strPath+"run038.root").c_str(),"hdat038","MuOnia, L=19.6/pb",3));

    //vecHistos.push_back(strctHisto::strctHisto((strPath+"run039.root").c_str(),"hdat039","MuOnia, L=19.6/pb",2));
    //vecHistos.push_back(strctHisto::strctHisto((strPath+"run043.root").c_str(),"hdat043","MuOnia, L=36.6/pb",3));
    //vecHistos.push_back(strctHisto::strctHisto((strPath+"run053.root").c_str(),"hdat053","MuOnia Trackrefit, L=39.9/pb",4));

    //vecHistos.push_back(strctHisto::strctHisto((strPath+"run053.root").c_str(),"hdat053","MuOnia L=41/pb refit",3));
    //vecHistos.push_back(strctHisto::strctHisto((strPath+"run074.root").c_str(),"hdat074","MuOnia L=41/pb cand 6",1));
    //vecHistos.push_back(strctHisto::strctHisto((strPath+"run075.root").c_str(),"hdat075","MuOnia L=41/pb cand 7",2));

//    vecHistos.push_back(strctHisto::strctHisto((strPath+"run078.root").c_str(),"hdat078","MuOnia L=41/pb ReReco cand 6",4));
//    vecHistos.push_back(strctHisto::strctHisto((strPath+"run080.root").c_str(),"hdat080","MCsig cand6match9",6));
//    vecHistos.push_back(strctHisto::strctHisto((strPath+"run081.root").c_str(),"hdat081","MCsig cand6",7));

//    vecHistos.push_back(strctHisto::strctHisto((strPath+"run082.root").c_str(),"hdat082","data cand6",2));
//    vecHistos.push_back(strctHisto::strctHisto((strPath+"run083.root").c_str(),"hdat083","data cand6match9",3));
//    vecHistos.push_back(strctHisto::strctHisto((strPath+"run086.root").c_str(),"hdat086","data cand7",4));

    //vecHistos.push_back(strctHisto::strctHisto((strPath+"run087.root").c_str(),"hdat087","data cand7match9",6));
    vecHistos.push_back(strctHisto::strctHisto((strPath+"run089.root").c_str(),"hdat089","MCsig cand7match9",7));

    //vecHistos.push_back(strctHisto::strctHisto("/scratch/frmeier/run016.root","hdat7","Data candidate 7",6));
    //vecHistos.push_back(strctHisto::strctHisto((strPath+"run041.root").c_str(),"hdat041","MC 100k evt 1st run",2));
    //vecHistos.push_back(strctHisto::strctHisto((strPath+"run042.root").c_str(),"hdat042","MC 100k evt 2nd run",3));


    //vecHistos.push_back(strctHisto::strctHisto((strPath+"run041.root").c_str(),"hdat041","MC 100k evt 1st run",2));
    //vecHistos.push_back(strctHisto::strctHisto((strPath+"run042.root").c_str(),"hdat042","MC 100k evt 2nd run",3));
    //vecHistos.push_back(strctHisto::strctHisto((strPath+"jpsiTest_signalMC_43_1_zlB.lambdaReader.default605122.cut.root").c_str(),"test","test",1));

    // check if the files were opened correctly
    {
	strctHisto_type::iterator it = vecHistos.begin();
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

    // Canvas
    c1 = new TCanvas;
    c1->SetWindowSize(1200,600);
    //c1 = new TCanvas("c1","Superimposed histos",1200,600);
    TPad* pad1= (TPad*)c1->cd();
    pad1->cd();
    pad1->SetTicks(0,0);
    pad1->SetRightMargin(0.20);
    pad1->SetTopMargin(0.10);

    // check no of histos and complain if problems may arise
    int nHistos = vecHistos.size();
    if (nHistos > 5) {
	cout << "Warning: You request to superimpose more than 5 histos." << endl
	   << "I'm running out of space for the statboxes..." << endl;
    }

    // plotting
    std::string strHist("("+toString(nBins)+","+toString(mymin)+","+toString(mymax)+")");
    cout << strHist << endl;
    double glMax = 0.; // Max entry of all histos
    for (strctHisto_type::iterator it = vecHistos.begin(); it!= vecHistos.end(); it++)
    {
	// draw, but the first one needs no "sames"
	cout << "Plotting histo " << (*it).name << endl;
	if (it == vecHistos.begin())
	{
	    (*it).tree->Draw((strPlot+">>"+(*it).name+strHist).c_str(), strCut.c_str());
	}
	else
	{
	    (*it).tree->Draw((strPlot+">>"+(*it).name+strHist).c_str(), strCut.c_str(),"sames");
	}
	(*it).histo = (TH1F*)gDirectory->GetList()->FindObject((*it).name.c_str());
	(*it).histo->SetLineColor((*it).color);
	// normalize if required
	if(doNorm) (*it).histo->Scale(1./(*it).histo->GetSumOfWeights());
	// collect max entry
	if (it == vecHistos.begin())
	{
	    glMax = (*it).histo->GetMaximum();
	}
	else
	{
	    glMax = (*it).histo->GetMaximum() > glMax ? (*it).histo->GetMaximum() : glMax;
	}
    }
    // scale max and set it below
    glMax *= 1.1; 
    for (strctHisto_type::iterator it = vecHistos.begin(); it!= vecHistos.end(); it++)
    {
	(*it).histo->SetMaximum(glMax);
    }
    
    // Statboxen
    pad1->Update();
    double deltay=0, top_corner = 0.9;
    double dHeightStatbox = .5/(double)nHistos;
    TPaveStats *st1;
    for (strctHisto_type::iterator it = vecHistos.begin(); it!= vecHistos.end(); it++)
    {
	st1 = (TPaveStats*)(*it).histo->GetListOfFunctions()->FindObject("stats");
        st1->SetOptStat(112211);
	st1->SetY1NDC(top_corner-deltay-dHeightStatbox);
        st1->SetY2NDC(top_corner-deltay);
	st1->SetX1NDC(0.80); st1->SetX2NDC(.995);
        st1->SetTextColor((*it).color);
	//st1->SetStatFormat("2.3e");
	//st1->SetFitFormat("2.3e");
	deltay+=dHeightStatbox;
    }

    // Legende einfügen
    TLegend *legend = new TLegend(0.80,top_corner-deltay-dHeightStatbox,0.995,top_corner-deltay);
    legend->SetFillStyle(1000);
    legend->SetBorderSize(1.);
    legend->SetTextSize(0.03);
    legend->SetFillColor(0);
    //legend->SetBorderSize(0.);
    for (strctHisto_type::iterator it = vecHistos.begin(); it!= vecHistos.end(); it++)
    {
	legend->AddEntry((*it).histo,(*it).title.c_str(),"pl");
    }
    legend->Draw();

    c1->SaveAs("c1.pdf");
}

