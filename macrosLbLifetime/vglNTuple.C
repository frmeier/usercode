#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"

TFile *f1, *f2;
TTree *t1, *t2, *t;
void vglNTuple()
{
    gSystem->Load("../../../CMSSW_3_8_4_patch2/src/AnalysisDataFormats/HeavyFlavorObjects/lib/libAna00.so");
    //f1 = TFile::Open("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/frmeier/lambda/l0044/jpsiTest_signalMC_1_1_Clt.root");
    //f2 = TFile::Open("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/frmeier/lambda/l0045/jpsiTest_signalMC_1_1_Yes.root");
    f1 = TFile::Open("/scratch/frmeier/jpsiTest_signalMC_1_1_Clt.root"); // l0044
    //f2 = TFile::Open("/scratch/frmeier/jpsiTest_signalMC_1_1_Yes.root"); // l0045
    t1 = (TTree*)f1->Get("T1");
    //t2 = (TTree*)f2->Get("T1");
    //t->AddFriend("T1",f1);
    //t1->AddFriend("T2=T1",f2);
    t1->AddFriend("T2=T1","/scratch/frmeier/jpsiTest_signalMC_1_1_Yes.root");
}

