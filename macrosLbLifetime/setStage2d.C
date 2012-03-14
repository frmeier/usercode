#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TPaletteAxis.h"
#include "TAxis.h"
#include "TMath.h"
#include "setTDRStyle_modified.C"
#include <string>
#include <sstream>
#include <iostream>
#include <vector>

using std::cout;
using std::endl;

TCanvas *c1;

struct strctHisto2d{
    strctHisto2d(std::string pFilename, std::string pName, std::string pTitle, int pColor)
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
    TH2F* histo;
    bool valid;
};

typedef std::vector<strctHisto2d> strctHisto2d_type;

template <typename T>
std::string toString(T i)
{
    ostringstream oss;
    oss << i;
    return oss.str();
}

void setStage2d(std::string strPlot, std::string strCut, int nBinsX, double myminX, double mymaxX, int nBinsY, double myminY, double mymaxY) {
    setTDRStyle();
    gStyle->SetPalette(1);

    gStyle->SetOptStat(111111);

    // Set up vector with histos
    strctHisto2d_type vecHistos;
    std::string strPath = "/scratch/frmeier/";
    //std::string strPath = "/home/frank/psi/lambda/data/";
    //std::string strPath = "/media/DR-1/";
    //std::string strPath = "dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/frmeier/readerRuns/mergedTrees/";
    //vecHistos.push_back(strctHisto2d::strctHisto2d((strPath+"run009.root").c_str(),"hsig6","Signal MC Cand 6",1));
    //vecHistos.push_back(strctHisto2d::strctHisto2d((strPath+"run012.root").c_str(),"hsig7","Signal MC Cand 7",2));
    //vecHistos.push_back(strctHisto2d::strctHisto2d((strPath+"run015.root").c_str(),"hsig8","Signal MC Cand 8",3));
    //vecHistos.push_back(strctHisto2d::strctHisto2d((strPath+"run018.root").c_str(),"hsig9","Signal MC Cand 9",4));
    //vecHistos.push_back(strctHisto2d::strctHisto2d((strPath+"run018.root").c_str(),"hsig","Signal MC Cand ",4));
    //vecHistos.push_back(strctHisto2d::strctHisto2d("/scratch/frmeier/run011part.root","hbgr6","Background MC",2));
    //vecHistos.push_back(strctHisto2d::strctHisto2d("/scratch/frmeier/run010.root","hdat6","Data",4));
    //vecHistos.push_back(strctHisto2d::strctHisto2d("/scratch/frmeier/run016.root","hdat7","Data candidate 7",6));
    //vecHistos.push_back(strctHisto2d::strctHisto2d((strPath+"run033.root").c_str(),"hdat033","MuOnia, L=12.2/pb",1));
    //vecHistos.push_back(strctHisto2d::strctHisto2d((strPath+"run009part.root").c_str(),"hsig009p","MC signal",1));
    //vecHistos.push_back(strctHisto2d::strctHisto2d((strPath+"run033.root").c_str(),"hdat033","MuOnia, L=12.2/pb",1));
    //vecHistos.push_back(strctHisto2d::strctHisto2d((strPath+"run034.root").c_str(),"hdat034","MuOnia, L=14.3/pb",4));
    vecHistos.push_back(strctHisto2d::strctHisto2d((strPath+"run038.root").c_str(),"hdat038","MuOnia, L=19.6/pb",3));
    vecHistos.push_back(strctHisto2d::strctHisto2d((strPath+"run039.root").c_str(),"hdat039","MuOnia, L=19.6/pb",2));

    // check if the files were opened correctly
    {
	strctHisto2d_type::iterator it = vecHistos.begin();
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
    c1->SetWindowSize(1000,1200);
    //c1 = new TCanvas("c1"," histos",1200,600);

    // check no of histos and complain if problems may arise
    int nHistos = vecHistos.size();
    if (nHistos > 9) {
	cout << "Warning: You request to draw more than 9 histos." << endl
	   << "I'm running out of space..." << endl;
    }
    if (nHistos == 0)
    {
	cout << "No histograms in list. Ending." << endl;
	return;
    }

    // calculate an 'optimal' size of the canvas
    int cSizeY = TMath::Ceil(sqrt(nHistos));
    int cSizeX = TMath::Ceil((double)nHistos/(double)cSizeY);
    c1->Divide(cSizeX,cSizeY);

    // plotting
    //std::string strHist("("+toString(nBinsX)+","+toString(myminX)+","+toString(mymaxX)+"," +toString(nBinsY)+","+toString(myminY)+","+toString(mymaxY)+")");
    std::string strHist("("+toString(nBinsY)+","+toString(myminY)+","+toString(mymaxY)+"," +toString(nBinsX)+","+toString(myminX)+","+toString(mymaxX)+")");
    cout << strHist << endl;

    //double glMax = 0.; // Max entry of all histos
    int cntHisto=1;
    for (strctHisto2d_type::iterator it = vecHistos.begin(); it!= vecHistos.end(); it++)
    {
	c1->cd(cntHisto);

	(*it).tree->Draw((strPlot+">>"+(*it).name+strHist).c_str(), strCut.c_str(),"COLZ");
	(*it).histo = (TH2F*)gDirectory->GetList()->FindObject((*it).name.c_str());
	//(*it).histo->SetLineColor((*it).color);

	// Make pad a bit smaller to have more space for the palette
	TPad* pad1= (TPad*)c1->cd(cntHisto);
	pad1->cd();
	pad1->SetRightMargin(0.20);
	pad1->SetLeftMargin(0.20);
	pad1->SetTopMargin(0.10);
	pad1->Update();
	// Reposition and resize palette
	TPaletteAxis *pal;
	pal = (TPaletteAxis*)(*it).histo->GetListOfFunctions()->FindObject("palette");
	pal->SetX1NDC(0.83);
	pal->SetY1NDC(0.14);
	pal->SetX2NDC(0.88);
	pal->SetY2NDC(0.60);
	pal->SetLabelSize(.04);
	// Add a separate TPaveText for the title
	TPaveText *ptTitle = new TPaveText(0.25,0.90,0.80,0.95,"LNDC");
	ptTitle->AddText((*it).title.c_str());
	ptTitle->SetBorderSize(0);
	ptTitle->SetFillColor(0);
	ptTitle->SetTextSize(0.04);
	ptTitle->Draw();
	// Adjust axis
	std::string titleY = strPlot.substr(0,strPlot.find(":"));
	std::string titleX = strPlot.substr(strPlot.find(":")+1,strPlot.size()-strPlot.find(":")-1);
	pad1->Update();
	(*it).histo->GetXaxis()->SetNdivisions(502);
	(*it).histo->GetXaxis()->SetTitle(titleX.c_str());
	(*it).histo->GetYaxis()->SetNdivisions(508);
	(*it).histo->GetYaxis()->SetTitle(titleY.c_str());
	(*it).histo->GetYaxis()->SetTitleOffset(1.5);
	cntHisto++;
    }
    /*
    // Statboxen
    pad1->Update();
    double deltay=0, top_corner = 0.9;
    double dHeightStatbox = .5/(double)nHistos;
    TPaveStats *st1;
    for (strctHisto2d_type::iterator it = vecHistos.begin(); it!= vecHistos.end(); it++)
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
    */

    /*
    // Legende einfügen
    TLegend *legend = new TLegend(0.80,top_corner-deltay-dHeightStatbox,0.995,top_corner-deltay);
    legend->SetFillStyle(1000);
    legend->SetBorderSize(1.);
    legend->SetTextSize(0.03);
    legend->SetFillColor(0);
    //legend->SetBorderSize(0.);
    for (strctHisto2d_type::iterator it = vecHistos.begin(); it!= vecHistos.end(); it++)
    {
	legend->AddEntry((*it).histo,(*it).title.c_str(),"pl");
    }
    legend->Draw();
    */
    c1->SaveAs("c1.pdf");
}

