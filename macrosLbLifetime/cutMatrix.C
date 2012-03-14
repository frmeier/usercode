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

struct strctCutset{
    strctCutset(std::string pCutLeft, std::vector<std::string> pCutsRight)
    {
	cutLeft = pCutLeft;
	cutsRight = pCutsRight;
    }
    strctCutset(std::string pCutLeft, const std::string pCutsRight[], unsigned int n)
    {
	cutLeft = pCutLeft;
	cutsRight.insert(cutsRight.begin(),pCutsRight,pCutsRight+n);
    }
    std::string getCut(unsigned int cutNo)
    {
	if (cutNo < cutsRight.size())
	{
	    return "&&" + cutLeft + cutsRight[cutNo];
	}
	else
	{
	    // if an inexistent item is requested, give something that does not harm
	    return "";
	}
    }
    unsigned int size()
    {
	return cutsRight.size();
    }
    std::string cutLeft;
    std::vector<std::string> cutsRight;
};

typedef std::vector<strctHisto> strctHisto_type;

template <typename T>
std::string toString(T i)
{
    ostringstream oss;
    oss << i;
    return oss.str();
}

void plotHistos(TPad* pad1, strctHisto_type vecHistos, std::string strPlot, std::string strCut, std::string strTitle,
		const int nBins, const double mymin, const double mymax, const bool doNorm = false);

void cutMatrix(std::string strPlot, int nBins, double mymin, double mymax, bool doNorm = false) {
	// set style
    setTDRStyle();
    gStyle->SetOptStat(112211);
    // Canvas
    c1 = new TCanvas;
    c1->SetWindowSize(1200,600);

    // Set up vector with histos
    strctHisto_type vecHistos;
    std::string strPath = "/scratch/frmeier/";
    //std::string strPath = "/home/frank/psi/lambda/daten/";
    //std::string strPath = "/media/DR-1/";

    vecHistos.push_back(strctHisto::strctHisto((strPath+"run019.root").c_str(),"hdat019","Data old",4));
    vecHistos.push_back(strctHisto::strctHisto((strPath+"run033.root").c_str(),"hdat033","Data new",1));

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
		std::cout << "The file for the histo " << it->name
		    << " " << it->title
		    << " could not be opened - histo dropped from list."
		    << std::endl;
		it = vecHistos.erase(it);
	    }
	}
    }

	// prepare the cuts and plot
    std::string constCut("(id1m&4)==4&&(id2m&4)==4&mjp>3.03&&mjp<3.16&&ptpr>ptpi&&alphal0<0.002&&ml0<1.15&&ptl0>4.&&d3lb/d3Elb>1.");

    std::vector<strctCutset> cutset;

    /*
    {
	std::string cuts[] = {"<0.002","<0.001","<0.0005","<0.004","<0.008"};
	cutset.push_back(strctCutset::strctCutset("alphal0",cuts,5));
    }
    {
	std::string cuts[] = {"<1.15","<1.14","<1.13","<1.16","<1.17"};
	cutset.push_back(strctCutset::strctCutset("ml0",cuts,5));
    }
    {
	std::string cuts[] = {">4.0",">2.0",">1.0",">8.0"};
	cutset.push_back(strctCutset::strctCutset("ptl0",cuts,4));
    }
    */
    {
	std::string cuts[] = {">5.0",">2.0",">3.5",">7.5",">10.0"};
	cutset.push_back(strctCutset::strctCutset("d3l0",cuts,5));
    }
    {
	std::string cuts[] = {">20.0",">10.0",">30",">40"};
	cutset.push_back(strctCutset::strctCutset("d3l0/d3El0",cuts,4));
    }

    // determine the max number of cuts
    unsigned int maxCutNo(0), noOfHistos(0);
    for(unsigned int i=0; i!=cutset.size(); i++)
	{
		noOfHistos+= cutset[i].size();
		if(cutset[i].size() > maxCutNo)
		{
			maxCutNo = cutset[i].size();
		}
	}
	std::cout << "We have " << noOfHistos << " histograms to draw" << std::endl;
	std::cout << "There are " << cutset.size() << " cutsets" << std::endl;
	std::cout << "and the max number of cuts per set is " << maxCutNo << std::endl;

	// now divide the canvas
	c1->Divide(maxCutNo, cutset.size());

    // now loop over all cuts and plot
    for(unsigned int i=0; i!=cutset.size(); i++)
    {
	std::cout << i << " - " << cutset[i].cutLeft << std::endl;
	for(unsigned int j=0; j!=maxCutNo; j++)
	{
	    std::string curCut(constCut);
	    std::string curTitle("");
	    bool cutsetOk(true);
	    for(unsigned int k=0; k!=cutset.size(); k++)
	    {
		if(i!=k)
		{
		    curCut+=cutset[k].getCut(0);
		}
		else
		{
		    if(j<cutset[k].size())
		    {
			curCut+=cutset[k].getCut(j);
			curTitle=cutset[k].getCut(j);
		    }
		    else
		    {
			cutsetOk = false;
		    }
		}
	    }
	    if(cutsetOk) // we can draw
	    {
			std::cout << curCut << std::endl;
			const int padNo = 1+i*maxCutNo+j;
			TPad* pad1= (TPad*)c1->cd(padNo);
			std::cout << pad1 << std::endl;
			plotHistos(pad1, vecHistos, strPlot, curCut, curTitle, nBins, mymin, mymax, doNorm);
	    }
	}
	std::cout << "-----" << std::endl;
    }

    return;


    c1->SaveAs("c1.pdf");
}

void plotHistos(TPad* pad1, strctHisto_type vecHistos, std::string strPlot, std::string strCut, std::string strTitle,
		const int nBins, const double mymin, const double mymax, const bool doNorm)
{	
    pad1->cd();
    const int curPadNo = pad1->GetNumber();
    pad1->SetTicks(0,0);
    pad1->SetRightMargin(0.20);
    pad1->SetTopMargin(0.10);

    // check no of histos and complain if problems may arise
    int nHistos = vecHistos.size();
    if (nHistos > 5) {
	std::cout << "Warning: You request to superimpose more than 5 histos." << std::endl
	   << "I'm running out of space for the statboxes..." << std::endl;
    }

    // plotting
    std::string strHist("("+toString(nBins)+","+toString(mymin)+","+toString(mymax)+")");
    std::cout << strHist << std::endl;
    double glMax = 0.; // Max entry of all histos
    for (strctHisto_type::iterator it = vecHistos.begin(); it!= vecHistos.end(); it++)
    {
	std::string curHistoName = (*it).name + toString(curPadNo);
	// draw, but the first one needs no "sames"
	if (it == vecHistos.begin())
	{
	    (*it).tree->Draw((strPlot+">>"+curHistoName+strHist).c_str(), strCut.c_str());
	}
	else
	{
	    (*it).tree->Draw((strPlot+">>"+curHistoName+strHist).c_str(), strCut.c_str(),"sames");
	}
	(*it).histo = (TH1F*)gDirectory->GetList()->FindObject(curHistoName.c_str());
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
	(*it).histo->SetTitle(strTitle.c_str());
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
	pad1->Update();
}

