#include "TROOT.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TAxis.h"
#include "utils.h"
#include <string>
#include <vector>

using std::string;

void doGraphError(std::vector<double> const &x, std::vector<double> const &y, std::vector<double> const &y_err,
	string title, string titleX, string unitX, string titleY, string unitY, double minY = -1111, double maxY = -1111)
{
    gStyle->SetPalette(1);
    if ((x.size() != y.size()) || (x.size() != y_err.size()) )
    {
	cout << "doGraphError: vectors do not match in size -- abort" << endl;
	return;
    }
    TGraphErrors *tge = new TGraphErrors(x.size(), &x[0], &y[0], 0, &y_err[0]);
    tge->Draw("AP");
    gPad->SetLeftMargin(.15);
    gPad->SetRightMargin(.2);
    gPad->SetTopMargin(.1);
    gPad->SetBottomMargin(.12);
    tge->SetTitle(title.c_str());
    tge->GetXaxis()->SetTitle(unitX.size()>0 ? (titleX+" / "+unitX).c_str() : titleX.c_str());
    //tge->GetXaxis()->SetNdivisions(505);
    tge->GetYaxis()->SetTitle(unitY.size()>0 ? (titleY+" / "+unitY).c_str() : titleY.c_str());
    tge->SetMinimum(minY);
    tge->SetMaximum(maxY);
    return; 
}

