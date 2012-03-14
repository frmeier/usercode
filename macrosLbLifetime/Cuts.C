#include <vector>
#include <string>
#include "Cuts.h"
#include "Cut.h"

// A collection of cuts

void Cuts::select1Cut(std::string name)
{
    if (name == "nocuts") { selectCutTrivial(); return; };
    if (name == "101208") { selectCut101208(); return; }
    // analysis cuts
    if (name == "an01") { selectCut101208(); cout << "WARNING! Cut an01 is deprecated. Use a newer one instead." << endl; return; }
    if (name == "an02") { selectCut110131(); cout << "WARNING! Cut an02 is deprecated. Use a newer one instead." << endl; return; }
    if (name == "an03") { selectCut110203(); cout << "WARNING! Cut an03 is deprecated. Use a newer one instead." << endl; return; }
    if (name == "an04") { selectCut110214(); cout << "WARNING! Cut an04 is deprecated. Use a newer one instead." << endl; return; }
    if (name == "an05") { selectCut110303(); return; }
    if (name == "an05exp") { selectCut110303exp(); return; }
    if (name == "an06") { selectCut110413(); return; }
    if (name == "an03exp") { selectCut110203exp2(); cout << "WARNING! an03exp is an experimental cutset which may change." << endl; return; }
    if (name == "an04exp") { selectCut110222exp(); cout << "WARNING! an04exp is an experimental cutset which may change." << endl; return; }
    if (name == "lb07") { selectCut110525(); return; }
    if (name == "lb08") { selectCut120309(); return; }
    // analysis cuts for B0
    if (name == "B001") { selectCutB0_110628(); return; }
    if (name == "B001exp") { selectCutB0_110628exp(); return; }
    // acceptance cuts
    if (name == "acc01") { selectCutAcceptance110120(); cout << "WARNING! Cut acc01 is deprecated. Use acc03 instead." << endl; return; }
    if (name == "acc02") { selectCutAcceptance110131(); cout << "WARNING! Cut acc02 is deprecated." << endl; return; }
    if (name == "acc03") { selectCutAcceptance110210(); return; }
    if (name == "acc03_16") { selectCutAcceptance110210(); selectCutMuEta16(); return; }
    if (name == "acc03B0") { selectCutAcceptance110210B0(); return; }
    if (name == "acc04Lb") { selectCutAcceptance111215global50(); selectCutAcceptance111215_L0(); return; }
    if (name == "acc04B0") { selectCutAcceptance111215global50(); selectCutAcceptance111215_Ks(); return; }
    // trigger cuts
    if (name == "L1_0_01") { selectCutTrgL1_110120_0(); return; }
    if (name == "L1_3_01") { selectCutTrgL1_110120_3(); return; }
    if (name == "HLT_01") { selectCutTrgHLT_110120(); return; }
    if (name == "HLT_02") { selectCutTrgHLT_110131(); return; }
    if (name == "HLT_matched_01") { selectCutTrgMatch_110209(); return; }
    if (name == "HLT_matched_2011_01") { selectCutTrgMatch_2011_110413(); return; }
    if (name == "HLT_matched_2011_02") { selectCutTrgMatch_2011_110512(); return; }

    if (name == "HLT_matched") { selectCutTrgMatch(); return; }
    if (name == "HLT_jpsi") { selectCutTrgJpsi(); return; }
    if (name == "HLT_jpsiDispl") { selectCutTrgJpsiDispl(); return; }
    if (name == "HLT_jpsiBarrel") { selectCutTrgJpsiBarrel(); return; }
    // special purpose cuts
    if (name == "mlbWindow01") { selectCutMlbWindow_110215(); return; }
    if (name == "iso01") { selectCutIso120229(); return; }
    // bin cuts
    if (name == "ptlb_10_15") { selectCutPtLb_10_15(); return; }
    if (name == "ptlb_15_20") { selectCutPtLb_15_20(); return; }
    if (name == "ptlb_20_infty") { selectCutPtLb_20_infty(); return; }
    if (name == "ylb_00_05") { selectCutYlb_00_05(); return; }
    if (name == "ylb_05_12") { selectCutYlb_05_12(); return; }
    if (name == "ylb_12_22") { selectCutYlb_12_22(); return; }
    // MC study cuts
    if (name == "isMCmatch") { cs.addCut(new cutConst("isMCmatch==1")); parvec.push_back(0); return; }
    if (name == "isSig") { cs.addCut(new cutConst("isSig==1")); parvec.push_back(0); return; }
    // nothing found
    cout << "Error: Cut not found. Exiting." << endl;
    throw ("Error: Cut not found. Exiting.");
};

void Cuts::selectCutTrivial()
{
    cs.addCut(new cutSymWindow("mjp",3.095));   parvec.push_back(20); 
}

void Cuts::selectCut101208()
{
    // jp
    cs.addCut(new cutBitcheck("rid1m"));        parvec.push_back(4);
    cs.addCut(new cutBitcheck("rid2m"));        parvec.push_back(4);
    cs.addCut(new cutSymWindow("mjp",3.095));   parvec.push_back(0.200); // mjp 0.065
    cs.addCut(new cutBoundLower("prob1m"));     parvec.push_back(0.10);
    cs.addCut(new cutBoundLower("prob2m"));     parvec.push_back(0.10);
    cs.addCut(new cutBoundLower("ptjp"));       parvec.push_back(2.0);
    cs.addCut(new cutBoundLower("probjp"));     parvec.push_back(0.005);
    // l0
    cs.addCut(new cutSymWindow("ml0",1.115));   parvec.push_back(0.0140); // ml0 0.015
    cs.addCut(new cutBoundLower("probpr"));     parvec.push_back(0.020);
    cs.addCut(new cutBoundLower("probpi"));     parvec.push_back(0.020);
    cs.addCut(new cutGT2Var("rptpr","rptpi"));  parvec.push_back(0);
    cs.addCut(new cutBoundLower("ptl0"));       parvec.push_back(3.0);   // ptl0 4.0
    cs.addCut(new cutBoundLower("rptpr"));      parvec.push_back(1.0);
    cs.addCut(new cutBoundLower("rptpi"));      parvec.push_back(0.50);
    cs.addCut(new cutBoundLower("probl0"));     parvec.push_back(0.020);
    //cs.addCut(new cutBoundLower("ctl0"));       parvec.push_back(1.0);
    cs.addCut(new cutBoundUpper("alphal0"));    parvec.push_back(0.0060);   // alpha 0.002
    cs.addCut(new cutBoundLower("d3l0"));       parvec.push_back(1.0);   // d3l0 5.0
    cs.addCut(new cutBoundLower("d3l0/d3El0")); parvec.push_back(10);    // d3l0/d3El0 30.0
    // lb
    //cs.addCut(new cutBoundLower("d3lb/d3Elb")); parvec.push_back(2.05);  // d3lb/d3Elb 1.0
    //cs.addCut(new cutBoundLower("ctlb"));       parvec.push_back(0.0020);
    cs.addCut(new cutBoundLower("problb"));     parvec.push_back(0.0010);
    cs.addCut(new cutBoundUpper("alphalb"));    parvec.push_back(0.30);   // alpha 0.002
}

void Cuts::selectCut110131()
{
    // jp
    cs.addCut(new cutBitcheck("rid1m"));        parvec.push_back(4);
    cs.addCut(new cutBitcheck("rid2m"));        parvec.push_back(4);
    cs.addCut(new cutSymWindow("mjp",3.095));   parvec.push_back(0.200); // mjp 0.065
    cs.addCut(new cutBoundLower("prob1m"));     parvec.push_back(0.10);
    cs.addCut(new cutBoundLower("prob2m"));     parvec.push_back(0.10);
    cs.addCut(new cutBoundLower("ptjp"));       parvec.push_back(2.0);
    cs.addCut(new cutBoundLower("probjp"));     parvec.push_back(0.005);
    // l0
    cs.addCut(new cutSymWindow("ml0",1.115));   parvec.push_back(0.0140); // ml0 0.015
    cs.addCut(new cutSymWindowVeto("Kshypo",0.497614));   parvec.push_back(0.008);
    //cs.addCut(new cutRatioWindow("iprpr",1.));  parvec.push_back(2.5);
    //cs.addCut(new cutRatioWindow("iprpi",1.));  parvec.push_back(2.5);
    //cs.addCut(new cutRatioWindow("ipr1m",1.));  parvec.push_back(2.5);
    //cs.addCut(new cutRatioWindow("ipr2m",1.));  parvec.push_back(2.5);
    cs.addCut(new cutBoundLower("probpr"));     parvec.push_back(0.020);
    cs.addCut(new cutBoundLower("probpi"));     parvec.push_back(0.020);
    cs.addCut(new cutGT2Var("rptpr","rptpi"));  parvec.push_back(0);
    cs.addCut(new cutBoundLower("ptl0"));       parvec.push_back(3.0);   // ptl0 4.0
    cs.addCut(new cutBoundLower("rptpr"));      parvec.push_back(1.0);
    cs.addCut(new cutBoundLower("rptpi"));      parvec.push_back(0.40);
    cs.addCut(new cutBoundLower("probl0"));     parvec.push_back(0.020);
    cs.addCut(new cutBoundUpper("alphal0"));    parvec.push_back(0.0060);   // alpha 0.002
    cs.addCut(new cutBoundLower("d3l0"));       parvec.push_back(1.0);   // d3l0 5.0
    cs.addCut(new cutBoundLower("d3l0/d3El0")); parvec.push_back(10);    // d3l0/d3El0 30.0
    // lb
    cs.addCut(new cutBoundLower("problb"));     parvec.push_back(0.0010);
    cs.addCut(new cutBoundUpper("alphalb"));    parvec.push_back(0.10);   // alpha 0.002
}

void Cuts::selectCut110203()
{
    // jp
    cs.addCut(new cutBitcheck("rid1m"));        parvec.push_back(4);
    cs.addCut(new cutBitcheck("rid2m"));        parvec.push_back(4);
    cs.addCut(new cutSymWindow("mjp",3.095));   parvec.push_back(0.200); // mjp 0.065
    cs.addCut(new cutBoundLower("prob1m"));     parvec.push_back(0.10);
    cs.addCut(new cutBoundLower("prob2m"));     parvec.push_back(0.10);
    cs.addCut(new cutBoundLower("ptjp"));       parvec.push_back(2.0);
    cs.addCut(new cutBoundLower("probjp"));     parvec.push_back(0.005);
    // l0
    cs.addCut(new cutSymWindow("ml0",1.115));   parvec.push_back(0.0140); // ml0 0.015
    cs.addCut(new cutSymWindowVeto("Kshypo",0.497614));   parvec.push_back(0.008);
    cs.addCut(new cutBoundLower("probpr"));     parvec.push_back(0.020);
    cs.addCut(new cutBoundLower("probpi"));     parvec.push_back(0.020);
    cs.addCut(new cutGT2Var("rptpr","rptpi"));  parvec.push_back(0);
    cs.addCut(new cutBoundLower("ptl0"));       parvec.push_back(2.0);   // ptl0 4.0
    cs.addCut(new cutBoundLower("rptpr"));      parvec.push_back(1.0);
    cs.addCut(new cutBoundLower("rptpi"));      parvec.push_back(0.40);
    cs.addCut(new cutBoundLower("probl0"));     parvec.push_back(0.020);
    cs.addCut(new cutBoundUpper("alphal0"));    parvec.push_back(0.0060);   // alpha 0.002
    cs.addCut(new cutBoundLower("d3l0"));       parvec.push_back(1.0);   // d3l0 5.0
    cs.addCut(new cutBoundLower("d3l0/d3El0")); parvec.push_back(10);    // d3l0/d3El0 30.0
    // lb
    cs.addCut(new cutBoundLower("problb"));     parvec.push_back(0.0010);
    cs.addCut(new cutBoundUpper("alphalb"));    parvec.push_back(0.10);   // alpha 0.002
}

void Cuts::selectCut110214()
{
    // jp
    cs.addCut(new cutBitcheck("rid1m"));        parvec.push_back(16);
    cs.addCut(new cutBitcheck("rid2m"));        parvec.push_back(16);
    cs.addCut(new cutSymWindow("mjp",3.095));   parvec.push_back(0.200); // mjp 0.065
    cs.addCut(new cutBoundLower("prob1m"));     parvec.push_back(0.10);
    cs.addCut(new cutBoundLower("prob2m"));     parvec.push_back(0.10);
    cs.addCut(new cutBoundLower("ptjp"));       parvec.push_back(2.0);
    cs.addCut(new cutBoundLower("probjp"));     parvec.push_back(0.005);
    // l0
    cs.addCut(new cutSymWindow("ml0",1.115));   parvec.push_back(0.0140); // ml0 0.015
    cs.addCut(new cutSymWindowVeto("Kshypo",0.497614));   parvec.push_back(0.008);
    cs.addCut(new cutBoundLower("probpr"));     parvec.push_back(0.020);
    cs.addCut(new cutBoundLower("probpi"));     parvec.push_back(0.020);
    cs.addCut(new cutGT2Var("rptpr","rptpi"));  parvec.push_back(0);
    cs.addCut(new cutBoundLower("ptl0"));       parvec.push_back(2.0);   // ptl0 4.0
    cs.addCut(new cutBoundLower("rptpr"));      parvec.push_back(1.0);
    cs.addCut(new cutBoundLower("rptpi"));      parvec.push_back(0.40);
    cs.addCut(new cutBoundLower("probl0"));     parvec.push_back(0.020);
    cs.addCut(new cutBoundUpper("alphal0"));    parvec.push_back(0.0060);   // alpha 0.002
    cs.addCut(new cutBoundLower("d3l0"));       parvec.push_back(1.0);   // d3l0 5.0
    cs.addCut(new cutBoundLower("d3l0/d3El0")); parvec.push_back(10);    // d3l0/d3El0 30.0
    // lb
    cs.addCut(new cutBoundLower("problb"));     parvec.push_back(0.0010);
    cs.addCut(new cutBoundUpper("alphalb"));    parvec.push_back(0.10);   // alpha 0.002
}

void Cuts::selectCut110303()
{
    // new cuts on ptlb and |ylb| to match kinematic range of final state in analysis
    // jp
    cs.addCut(new cutBitcheck("rid1m"));        parvec.push_back(16);
    cs.addCut(new cutBitcheck("rid2m"));        parvec.push_back(16);
    cs.addCut(new cutSymWindow("mjp",3.095));   parvec.push_back(0.200); // mjp 0.065
    cs.addCut(new cutBoundLower("prob1m"));     parvec.push_back(0.10);
    cs.addCut(new cutBoundLower("prob2m"));     parvec.push_back(0.10);
    cs.addCut(new cutBoundLower("ptjp"));       parvec.push_back(2.0);
    cs.addCut(new cutBoundLower("probjp"));     parvec.push_back(0.005);
    // l0
    cs.addCut(new cutSymWindow("ml0",1.115));   parvec.push_back(0.0140); // ml0 0.015
    cs.addCut(new cutSymWindowVeto("Kshypo",0.497614));   parvec.push_back(0.008);
    cs.addCut(new cutBoundLower("probpr"));     parvec.push_back(0.020);
    cs.addCut(new cutBoundLower("probpi"));     parvec.push_back(0.020);
    cs.addCut(new cutGT2Var("rptpr","rptpi"));  parvec.push_back(0);
    cs.addCut(new cutBoundLower("ptl0"));       parvec.push_back(2.0);   // ptl0 4.0
    cs.addCut(new cutBoundLower("rptpr"));      parvec.push_back(1.0);
    cs.addCut(new cutBoundLower("rptpi"));      parvec.push_back(0.40);
    cs.addCut(new cutBoundLower("probl0"));     parvec.push_back(0.020);
    cs.addCut(new cutBoundUpper("alphal0"));    parvec.push_back(0.0060);   // alpha 0.002
    cs.addCut(new cutBoundLower("d3l0"));       parvec.push_back(1.0);   // d3l0 5.0
    cs.addCut(new cutBoundLower("d3l0/d3El0")); parvec.push_back(10);    // d3l0/d3El0 30.0
    // lb
    cs.addCut(new cutBoundLower("problb"));     parvec.push_back(0.0010);
    cs.addCut(new cutBoundUpper("alphalb"));    parvec.push_back(0.10);   // alpha 0.002
    cs.addCut(new cutBoundLower("ptlb"));       parvec.push_back(6.0);
    cs.addCut(new cutConst("TMath::Abs(ylb)<2.2")); parvec.push_back(0);
}

void Cuts::selectCut110303exp()
{
    // new cuts on ptlb and |ylb| to match kinematic range of final state in analysis
    // jp
    cs.addCut(new cutBitcheck("rid1m"));        parvec.push_back(16);
    cs.addCut(new cutBitcheck("rid2m"));        parvec.push_back(16);
    cs.addCut(new cutSymWindow("mjp",3.095));   parvec.push_back(0.10); // mjp 0.065
    cs.addCut(new cutBoundLower("prob1m"));     parvec.push_back(0.10);
    cs.addCut(new cutBoundLower("prob2m"));     parvec.push_back(0.10);
    cs.addCut(new cutBoundLower("ptjp"));       parvec.push_back(2.0);
    //cs.addCut(new cutBoundLower("probjp"));     parvec.push_back(0.005);
    // l0
    cs.addCut(new cutSymWindow("ml0",1.115));   parvec.push_back(0.0050); // ml0 0.015
    cs.addCut(new cutSymWindowVeto("Kshypo",0.497614));   parvec.push_back(0.010);
    cs.addCut(new cutBoundLower("probpr"));     parvec.push_back(0.020);
    cs.addCut(new cutBoundLower("probpi"));     parvec.push_back(0.020);
    cs.addCut(new cutGT2Var("rptpr","rptpi"));  parvec.push_back(0);
    cs.addCut(new cutBoundLower("ptl0"));       parvec.push_back(2.2);   // ptl0 4.0
    cs.addCut(new cutBoundLower("rptpr"));      parvec.push_back(1.0);
    cs.addCut(new cutBoundLower("rptpi"));      parvec.push_back(0.50);
    //cs.addCut(new cutBoundLower("probl0"));     parvec.push_back(0.020);
    //cs.addCut(new cutBoundUpper("alphal0"));    parvec.push_back(0.0100);   // alpha 0.002
    cs.addCut(new cutBoundUpper("ptgangDRl0"));    parvec.push_back(0.02);   // Extra gross fuer cutscan
    cs.addCut(new cutBoundLower("d3l0"));       parvec.push_back(1.0);   // d3l0 5.0
    cs.addCut(new cutBoundLower("d3l0/d3El0")); parvec.push_back(5);    // d3l0/d3El0 30.0
    // lb
    //cs.addCut(new cutBoundLower("problb"));     parvec.push_back(0.0010);
    //cs.addCut(new cutBoundUpper("alphalb"));    parvec.push_back(0.16);   // alpha 0.002
    cs.addCut(new cutBoundLower("ptlb"));       parvec.push_back(6.0);
    cs.addCut(new cutConst("TMath::Abs(ylb)<2.2")); parvec.push_back(0);
}

void Cuts::selectCut110413()
{
    // Cuts to improve S/B on 2011A data
    // jp
    cs.addCut(new cutBitcheck("rid1m"));        parvec.push_back(16);
    cs.addCut(new cutBitcheck("rid2m"));        parvec.push_back(16);
    cs.addCut(new cutSymWindow("mjp",3.095));   parvec.push_back(0.10); // mjp 0.065
    cs.addCut(new cutBoundLower("prob1m"));     parvec.push_back(0.10);
    cs.addCut(new cutBoundLower("prob2m"));     parvec.push_back(0.10);
    cs.addCut(new cutBoundLower("ptjp"));       parvec.push_back(2.0);
    cs.addCut(new cutBoundLower("probjp"));     parvec.push_back(0.005);
    // l0
    cs.addCut(new cutSymWindow("ml0",1.115));   parvec.push_back(0.0050); // ml0 0.015
    cs.addCut(new cutSymWindowVeto("Kshypo",0.497614));   parvec.push_back(0.010);
    cs.addCut(new cutBoundLower("probpr"));     parvec.push_back(0.020);
    cs.addCut(new cutBoundLower("probpi"));     parvec.push_back(0.020);
    cs.addCut(new cutGT2Var("rptpr","rptpi"));  parvec.push_back(0);
    cs.addCut(new cutBoundLower("ptl0"));       parvec.push_back(2.2);   // ptl0 4.0
    cs.addCut(new cutBoundLower("rptpr"));      parvec.push_back(1.0);
    cs.addCut(new cutBoundLower("rptpi"));      parvec.push_back(0.50);
    cs.addCut(new cutBoundLower("probl0"));     parvec.push_back(0.020);
    //cs.addCut(new cutBoundUpper("alphal0"));    parvec.push_back(0.0100);   // alpha 0.002
    cs.addCut(new cutBoundUpper("ptgangDRl0"));    parvec.push_back(0.02);   // Extra gross fuer cutscan
    cs.addCut(new cutBoundLower("d3l0"));       parvec.push_back(1.0);   // d3l0 5.0
    cs.addCut(new cutBoundLower("d3l0/d3El0")); parvec.push_back(5);    // d3l0/d3El0 30.0
    // lb
    cs.addCut(new cutBoundLower("problb"));     parvec.push_back(0.0010);
    cs.addCut(new cutBoundUpper("alphalb"));    parvec.push_back(0.16);   // alpha 0.002
    cs.addCut(new cutBoundLower("ptlb"));       parvec.push_back(6.0);
    cs.addCut(new cutConst("TMath::Abs(ylb)<2.2")); parvec.push_back(0);
}

void Cuts::selectCut110525()
{
    // Cuts reflecting latest changed to triggers
    // jp
    cs.addCut(new cutBitcheck("rid1m"));        parvec.push_back(16);
    cs.addCut(new cutBitcheck("rid2m"));        parvec.push_back(16);
    cs.addCut(new cutSymWindow("mjp",3.095));   parvec.push_back(0.10); // mjp 0.065
    cs.addCut(new cutBoundLower("prob1m"));     parvec.push_back(0.10);
    cs.addCut(new cutBoundLower("prob2m"));     parvec.push_back(0.10);
    cs.addCut(new cutBoundLower("rpt1m"));      parvec.push_back(3.0);
    cs.addCut(new cutBoundLower("rpt2m"));      parvec.push_back(3.0);
    cs.addCut(new cutBoundLower("ptjp"));       parvec.push_back(6.5);
    cs.addCut(new cutBoundLower("probjp"));     parvec.push_back(0.06);
    // l0
    cs.addCut(new cutSymWindow("ml0",1.115));   parvec.push_back(0.0050); // ml0 0.015
    cs.addCut(new cutSymWindowVeto("Kshypo",0.497614));   parvec.push_back(0.010);
    cs.addCut(new cutBoundLower("probpr"));     parvec.push_back(0.020);
    cs.addCut(new cutBoundLower("probpi"));     parvec.push_back(0.020);
    cs.addCut(new cutGT2Var("rptpr","rptpi"));  parvec.push_back(0);
    cs.addCut(new cutBoundLower("ptl0"));       parvec.push_back(2.2);   // ptl0 4.0
    cs.addCut(new cutBoundLower("rptpr"));      parvec.push_back(1.8);
    cs.addCut(new cutBoundLower("rptpi"));      parvec.push_back(0.30);
    cs.addCut(new cutBoundLower("probl0"));     parvec.push_back(0.06);
    //cs.addCut(new cutBoundUpper("alphal0"));    parvec.push_back(0.0100);   // alpha 0.002
    cs.addCut(new cutBoundUpper("ptgangDRl0"));    parvec.push_back(0.03);   // Extra gross fuer cutscan
    cs.addCut(new cutBoundLower("d3l0"));       parvec.push_back(1.0);   // d3l0 5.0
    cs.addCut(new cutBoundLower("d3l0/d3El0")); parvec.push_back(5);    // d3l0/d3El0 30.0
    // lb
    cs.addCut(new cutBoundLower("problb"));     parvec.push_back(0.06);
    //cs.addCut(new cutBoundUpper("alphalb"));    parvec.push_back(0.16);   // alpha 0.002
    cs.addCut(new cutBoundLower("ptlb"));       parvec.push_back(6.0);
    cs.addCut(new cutConst("TMath::Abs(ylb)<2.2")); parvec.push_back(0);
}

void Cuts::selectCut120309()
{
    // Cuts after reviewing 2011 data
    // jp
    cs.addCut(new cutBitcheck("rid1m"));        parvec.push_back(16);
    cs.addCut(new cutBitcheck("rid2m"));        parvec.push_back(16);
    cs.addCut(new cutSymWindow("mjp",3.096916));   parvec.push_back(0.10); // mjp 0.065
    cs.addCut(new cutBoundLower("prob1m"));     parvec.push_back(0.10);
    cs.addCut(new cutBoundLower("prob2m"));     parvec.push_back(0.10);
    cs.addCut(new cutBoundLower("rpt1m"));      parvec.push_back(3.0);
    cs.addCut(new cutBoundLower("rpt2m"));      parvec.push_back(3.0);
    cs.addCut(new cutBoundLower("ptjp"));       parvec.push_back(6.5);
    cs.addCut(new cutBoundLower("probjp"));     parvec.push_back(0.06);
    // l0
    cs.addCut(new cutSymWindow("ml0",1.115683));   parvec.push_back(0.0050); // ml0 0.015
    cs.addCut(new cutSymWindowVeto("Kshypo",0.497614));   parvec.push_back(0.010);
    cs.addCut(new cutBoundLower("probpr"));     parvec.push_back(0.020);
    cs.addCut(new cutBoundLower("probpi"));     parvec.push_back(0.020);
    cs.addCut(new cutGT2Var("rptpr","rptpi"));  parvec.push_back(0);
    cs.addCut(new cutBoundLower("ptl0"));       parvec.push_back(2.2);   // ptl0 4.0
    cs.addCut(new cutBoundLower("rptpr"));      parvec.push_back(1.8);
    cs.addCut(new cutBoundLower("rptpi"));      parvec.push_back(0.30);
    cs.addCut(new cutBoundLower("probl0"));     parvec.push_back(0.06);
    //cs.addCut(new cutBoundUpper("alphal0"));    parvec.push_back(0.0100);   // alpha 0.002
    cs.addCut(new cutBoundUpper("ptgangDRl0"));    parvec.push_back(0.01);   // Extra gross fuer cutscan
    cs.addCut(new cutBoundLower("d3l0"));       parvec.push_back(1.0);   // d3l0 5.0
    cs.addCut(new cutBoundLower("d3l0/d3El0")); parvec.push_back(5);    // d3l0/d3El0 30.0
    // lb
    cs.addCut(new cutBoundLower("problb"));     parvec.push_back(0.06);
    //cs.addCut(new cutBoundUpper("alphalb"));    parvec.push_back(0.16);   // alpha 0.002
    //cs.addCut(new cutBoundLower("ptlb"));       parvec.push_back(6.0); // 2011 triggers make ptlb>9 and ptlb should not be a cut variable for lifetiem analysis
    //cs.addCut(new cutConst("TMath::Abs(ylb)<2.2")); parvec.push_back(0); // Barrel trigeer renders this useless...
}

void Cuts::selectCut110203exp()
{
    // jp
    cs.addCut(new cutBitcheck("rid1m"));        parvec.push_back(4);
    cs.addCut(new cutBitcheck("rid2m"));        parvec.push_back(4);
    cs.addCut(new cutSymWindow("mjp",3.095));   parvec.push_back(0.200); // mjp 0.065
    cs.addCut(new cutBoundLower("prob1m"));     parvec.push_back(0.10);
    cs.addCut(new cutBoundLower("prob2m"));     parvec.push_back(0.10);
    cs.addCut(new cutBoundLower("ptjp"));       parvec.push_back(2.0);
    cs.addCut(new cutBoundLower("probjp"));     parvec.push_back(0.005);
    // l0
    cs.addCut(new cutSymWindow("ml0",1.115));   parvec.push_back(0.0140); // ml0 0.015
    cs.addCut(new cutSymWindowVeto("Kshypo",0.497614));   parvec.push_back(0.016);
    cs.addCut(new cutBoundLower("probpr"));     parvec.push_back(0.020);
    cs.addCut(new cutBoundLower("probpi"));     parvec.push_back(0.020);
    cs.addCut(new cutGT2Var("rptpr","rptpi"));  parvec.push_back(0);
    cs.addCut(new cutBoundLower("ptl0"));       parvec.push_back(1.8);   // ptl0 4.0
    cs.addCut(new cutBoundLower("rptpr"));      parvec.push_back(0.8);
    cs.addCut(new cutBoundLower("rptpi"));      parvec.push_back(0.30);
    cs.addCut(new cutBoundLower("probl0"));     parvec.push_back(0.020);
    cs.addCut(new cutBoundUpper("alphal0"));    parvec.push_back(0.0040);   // alpha 0.002
    cs.addCut(new cutBoundLower("d3l0"));       parvec.push_back(1.0);   // d3l0 5.0
    cs.addCut(new cutBoundLower("d3l0/d3El0")); parvec.push_back(10);    // d3l0/d3El0 30.0
    // lb
    cs.addCut(new cutBoundLower("problb"));     parvec.push_back(0.0010);
    cs.addCut(new cutBoundUpper("alphalb"));    parvec.push_back(0.10);   // alpha 0.002
}

void Cuts::selectCut110203exp2()
{
    // jp
    cs.addCut(new cutSymWindow("mjp",3.095));   parvec.push_back(0.200); // mjp 0.065
    cs.addCut(new cutBoundLower("prob1m"));     parvec.push_back(0.10);
    cs.addCut(new cutBoundLower("prob2m"));     parvec.push_back(0.10);
    cs.addCut(new cutBoundLower("ptjp"));       parvec.push_back(2.0);
    cs.addCut(new cutBoundLower("probjp"));     parvec.push_back(0.005);
    // l0
    cs.addCut(new cutSymWindow("ml0",1.115));   parvec.push_back(0.0140); // ml0 0.015
    cs.addCut(new cutSymWindowVeto("Kshypo",0.497614));   parvec.push_back(0.008);
    cs.addCut(new cutBoundLower("probpr"));     parvec.push_back(0.020);
    cs.addCut(new cutBoundLower("probpi"));     parvec.push_back(0.020);
    cs.addCut(new cutGT2Var("rptpr","rptpi"));  parvec.push_back(0);
    cs.addCut(new cutBoundLower("ptl0"));       parvec.push_back(2.0);   // ptl0 4.0
    cs.addCut(new cutBoundLower("rptpr"));      parvec.push_back(1.0);
    cs.addCut(new cutBoundLower("rptpi"));      parvec.push_back(0.40);
    cs.addCut(new cutBoundLower("probl0"));     parvec.push_back(0.020);
    cs.addCut(new cutBoundUpper("alphal0"));    parvec.push_back(0.0060);   // alpha 0.002
    cs.addCut(new cutBoundLower("d3l0"));       parvec.push_back(1.0);   // d3l0 5.0
    cs.addCut(new cutBoundLower("d3l0/d3El0")); parvec.push_back(10);    // d3l0/d3El0 30.0
    // lb
    cs.addCut(new cutBoundLower("problb"));     parvec.push_back(0.0010);
    cs.addCut(new cutBoundUpper("alphalb"));    parvec.push_back(0.10);   // alpha 0.002
}

void Cuts::selectCut110222exp()
{
    // jp
    cs.addCut(new cutBitcheck("rid1m"));        parvec.push_back(16);
    cs.addCut(new cutBitcheck("rid2m"));        parvec.push_back(16);
    cs.addCut(new cutSymWindow("mjp",3.095));   parvec.push_back(0.200); // mjp 0.065
    cs.addCut(new cutBoundLower("prob1m"));     parvec.push_back(0.10);
    cs.addCut(new cutBoundLower("prob2m"));     parvec.push_back(0.10);
    cs.addCut(new cutBoundLower("ptjp"));       parvec.push_back(2.0);
    cs.addCut(new cutBoundLower("probjp"));     parvec.push_back(0.005);
    // l0
    cs.addCut(new cutSymWindow("ml0",1.115));   parvec.push_back(0.0140); // ml0 0.015
    cs.addCut(new cutSymWindowVeto("Kshypo",0.497614));   parvec.push_back(0.008);
    cs.addCut(new cutBoundLower("probpr"));     parvec.push_back(0.020);
    cs.addCut(new cutBoundLower("probpi"));     parvec.push_back(0.020);
    cs.addCut(new cutGT2Var("rptpr","rptpi"));  parvec.push_back(0);
    cs.addCut(new cutBoundLower("ptl0"));       parvec.push_back(2.0);   // ptl0 4.0
    cs.addCut(new cutBoundLower("rptpr"));      parvec.push_back(1.0);
    cs.addCut(new cutBoundLower("rptpi"));      parvec.push_back(0.40);
    cs.addCut(new cutBoundLower("probl0"));     parvec.push_back(0.020);
    //cs.addCut(new cutBoundUpper("alphal0"));    parvec.push_back(3.142);   // Extra gross fuer cutscan
    //cs.addCut(new cutBoundUpper("ptgangDRl0"));    parvec.push_back(12);   // Extra gross fuer cutscan
    cs.addCut(new cutBoundUpper("alphal0"));    parvec.push_back(0.0060);   // alpha 0.002
    cs.addCut(new cutBoundLower("d3l0"));       parvec.push_back(1.0);   // d3l0 5.0
    cs.addCut(new cutBoundLower("d3l0/d3El0")); parvec.push_back(10);    // d3l0/d3El0 30.0
    // lb
    cs.addCut(new cutBoundLower("problb"));     parvec.push_back(0.0010);
    //cs.addCut(new cutBoundUpper("alphalb"));    parvec.push_back(0.10);   // alpha 0.002
    cs.addCut(new cutBoundUpper("alphalb"));    parvec.push_back(3.141);   // Extra gross
    cs.addCut(new cutBoundUpper("ptgangDRlb"));    parvec.push_back(12);   // Extra gross fuer cutscan
}

// cuts for B0
void Cuts::selectCutB0_110628()
{
    // Cuts reflecting latest changed to triggers
    // jp
    cs.addCut(new cutBitcheck("rid1m"));        parvec.push_back(16);
    cs.addCut(new cutBitcheck("rid2m"));        parvec.push_back(16);
    cs.addCut(new cutSymWindow("mjp",3.095));   parvec.push_back(0.10); // mjp 0.065
    cs.addCut(new cutBoundLower("prob1m"));     parvec.push_back(0.10);
    cs.addCut(new cutBoundLower("prob2m"));     parvec.push_back(0.10);
    cs.addCut(new cutBoundLower("rpt1m"));      parvec.push_back(3.0);
    cs.addCut(new cutBoundLower("rpt2m"));      parvec.push_back(3.0);
    cs.addCut(new cutBoundLower("ptjp"));       parvec.push_back(6.5);
    cs.addCut(new cutBoundLower("probjp"));     parvec.push_back(0.06);
    // Ks
    cs.addCut(new cutSymWindow("mKs",0.4976));   parvec.push_back(0.0050); // ml0 0.015
    //cs.addCut(new cutSymWindowVeto("Kshypo",0.497614));   parvec.push_back(0.010);
    cs.addCut(new cutBoundLower("probpi1"));     parvec.push_back(0.020);
    cs.addCut(new cutBoundLower("probpi2"));     parvec.push_back(0.020);
    cs.addCut(new cutBoundLower("ptKs"));       parvec.push_back(2.2);   // ptl0 4.0
    cs.addCut(new cutBoundLower("rptpi1"));      parvec.push_back(0.50);
    cs.addCut(new cutBoundLower("rptpi2"));      parvec.push_back(0.50);
    cs.addCut(new cutBoundLower("probKs"));     parvec.push_back(0.06);
    //cs.addCut(new cutBoundUpper("alphal0"));    parvec.push_back(0.0100);   // alpha 0.002
    cs.addCut(new cutBoundUpper("ptgangDRKs"));    parvec.push_back(0.03);   // Extra gross fuer cutscan
    cs.addCut(new cutBoundLower("d3Ks"));       parvec.push_back(1.0);   // d3l0 5.0
    cs.addCut(new cutBoundLower("d3Ks/d3EKs")); parvec.push_back(5);    // d3l0/d3El0 30.0
    // lb
    cs.addCut(new cutBoundLower("probB0"));     parvec.push_back(0.06);
    //cs.addCut(new cutBoundUpper("alphalb"));    parvec.push_back(0.16);   // alpha 0.002
    cs.addCut(new cutBoundLower("ptB0"));       parvec.push_back(6.0);
    cs.addCut(new cutConst("TMath::Abs(yB0)<2.2")); parvec.push_back(0);
}

void Cuts::selectCutB0_110628exp()
{
    // Cuts reflecting latest changed to triggers
    // jp
    cs.addCut(new cutBitcheck("rid1m"));        parvec.push_back(16);
    cs.addCut(new cutBitcheck("rid2m"));        parvec.push_back(16);
    cs.addCut(new cutSymWindow("mjp",3.095));   parvec.push_back(0.10); // mjp 0.065
    cs.addCut(new cutBoundLower("prob1m"));     parvec.push_back(0.10);
    cs.addCut(new cutBoundLower("prob2m"));     parvec.push_back(0.10);
    cs.addCut(new cutBoundLower("rpt1m"));      parvec.push_back(3.0);
    cs.addCut(new cutBoundLower("rpt2m"));      parvec.push_back(3.0);
    cs.addCut(new cutBoundLower("ptjp"));       parvec.push_back(6.5);
    cs.addCut(new cutBoundLower("probjp"));     parvec.push_back(0.06);
    // Ks
    cs.addCut(new cutSymWindow("mKs",0.4976));   parvec.push_back(0.0050); // ml0 0.015
    //cs.addCut(new cutSymWindowVeto("Kshypo",0.497614));   parvec.push_back(0.010);
    cs.addCut(new cutBoundLower("probpi1"));     parvec.push_back(0.020);
    cs.addCut(new cutBoundLower("probpi2"));     parvec.push_back(0.020);
    cs.addCut(new cutBoundLower("ptKs"));       parvec.push_back(2.2);   // ptl0 4.0
    cs.addCut(new cutBoundLower("rptpi1"));      parvec.push_back(0.50);
    cs.addCut(new cutBoundLower("rptpi2"));      parvec.push_back(0.50);
    cs.addCut(new cutBoundLower("probKs"));     parvec.push_back(0.06);
    //cs.addCut(new cutBoundUpper("alphal0"));    parvec.push_back(0.0100);   // alpha 0.002
    cs.addCut(new cutBoundUpper("ptgangDRKs"));    parvec.push_back(0.0005);   // Extra gross fuer cutscan
    cs.addCut(new cutBoundLower("d3Ks"));       parvec.push_back(1.0);   // d3l0 5.0
    cs.addCut(new cutBoundLower("d3Ks/d3EKs")); parvec.push_back(5);    // d3l0/d3El0 30.0
    // lb
    cs.addCut(new cutBoundLower("probB0"));     parvec.push_back(0.06);
    //cs.addCut(new cutBoundUpper("alphalb"));    parvec.push_back(0.16);   // alpha 0.002
    cs.addCut(new cutBoundLower("ptB0"));       parvec.push_back(6.0);
    cs.addCut(new cutConst("TMath::Abs(yB0)<2.2")); parvec.push_back(0);
}


// this cut should never be modified by any studies, so it is a const cut
void Cuts::selectCutAcceptance110120()
{	// now deprecated, does only check on one muon, not on both!
    std::string accstr1 = "(TMath::Abs(Seta1m)<1.3&&Spt1m>3.3)";
    std::string accstr2 = "(TMath::Abs(Seta1m)>=1.3&&TMath::Abs(Seta1m)<2.2&&Sp1m>2.9)";
    std::string accstr3 = "(TMath::Abs(Seta1m)>=2.2&&TMath::Abs(Seta1m)<2.4&&Spt1m>0.8)";
    cs.addCut(new cutConst("(" + accstr1 + "||" + accstr2 +  "||" + accstr3 + ")"));
    parvec.push_back(0);
    cs.addCut(new cutConst("(vrl0>1&&vrl0<35&&vzl0<100)"));
    parvec.push_back(0);
}

void Cuts::selectCutAcceptance110131()
{
    std::string accstr1 = "(TMath::Abs(Seta1m)<1.3&&Spt1m>3.3)";
    std::string accstr2 = "(TMath::Abs(Seta1m)>=1.3&&TMath::Abs(Seta1m)<2.2&&Sp1m>2.9)";
    std::string accstr3 = "(TMath::Abs(Seta1m)>=2.2&&TMath::Abs(Seta1m)<2.4&&Spt1m>0.8)";
    cs.addCut(new cutConst("(" + accstr1 + "||" + accstr2 +  "||" + accstr3 + ")"));
    parvec.push_back(0);
    // change to above: vrl0 = sin(theta) for theta=11.5 deg (eta=2.4)
    cs.addCut(new cutConst("(vrl0>.18&&vrl0<35&&vzl0<100)"));
    parvec.push_back(0);
}

void Cuts::selectCutAcceptance110210()
{
    std::string accstr1m_1 = "(TMath::Abs(Seta1m)<1.3&&Spt1m>3.3)";
    std::string accstr1m_2 = "(TMath::Abs(Seta1m)>=1.3&&TMath::Abs(Seta1m)<2.2&&Sp1m>2.9)";
    std::string accstr1m_3 = "(TMath::Abs(Seta1m)>=2.2&&TMath::Abs(Seta1m)<2.4&&Spt1m>0.8)";
    cs.addCut(new cutConst("(" + accstr1m_1 + "||" + accstr1m_2 +  "||" + accstr1m_3 + ")"));
    parvec.push_back(0);
    std::string accstr2m_1 = "(TMath::Abs(Seta2m)<1.3&&Spt2m>3.3)";
    std::string accstr2m_2 = "(TMath::Abs(Seta2m)>=1.3&&TMath::Abs(Seta2m)<2.2&&Sp2m>2.9)";
    std::string accstr2m_3 = "(TMath::Abs(Seta2m)>=2.2&&TMath::Abs(Seta2m)<2.4&&Spt2m>0.8)";
    cs.addCut(new cutConst("(" + accstr2m_1 + "||" + accstr2m_2 +  "||" + accstr2m_3 + ")"));
    parvec.push_back(0);
    cs.addCut(new cutConst("(vrl0>1&&vrl0<35&&vzl0<100)"));
    parvec.push_back(0);
}

void Cuts::selectCutAcceptance110210B0()
{
    std::string accstr1m_1 = "(TMath::Abs(Seta1m)<1.3&&Spt1m>3.3)";
    std::string accstr1m_2 = "(TMath::Abs(Seta1m)>=1.3&&TMath::Abs(Seta1m)<2.2&&Sp1m>2.9)";
    std::string accstr1m_3 = "(TMath::Abs(Seta1m)>=2.2&&TMath::Abs(Seta1m)<2.4&&Spt1m>0.8)";
    cs.addCut(new cutConst("(" + accstr1m_1 + "||" + accstr1m_2 +  "||" + accstr1m_3 + ")"));
    parvec.push_back(0);
    std::string accstr2m_1 = "(TMath::Abs(Seta2m)<1.3&&Spt2m>3.3)";
    std::string accstr2m_2 = "(TMath::Abs(Seta2m)>=1.3&&TMath::Abs(Seta2m)<2.2&&Sp2m>2.9)";
    std::string accstr2m_3 = "(TMath::Abs(Seta2m)>=2.2&&TMath::Abs(Seta2m)<2.4&&Spt2m>0.8)";
    cs.addCut(new cutConst("(" + accstr2m_1 + "||" + accstr2m_2 +  "||" + accstr2m_3 + ")"));
    parvec.push_back(0);
    cs.addCut(new cutConst("(vrKs>1&&vrKs<35&&vzKs<100)"));
    parvec.push_back(0);
}

void Cuts::selectCutAcceptance111215tracker50() // muons based on AN-11-417, see Journal #13 p66 for calculation of cuts
{
    const std::string accstr1m_1 = "(TMath::Abs(Seta1m)<1.2&&Spt1m>3.5)";
    const std::string accstr1m_2 = "(TMath::Abs(Seta1m)>=1.2&&TMath::Abs(Seta1m)<1.6&&Sp1m>8-3.75*TMath::Abs(Seta1m))";
    const std::string accstr1m_3 = "(TMath::Abs(Seta1m)>=1.6&&TMath::Abs(Seta1m)<2.4&&Spt1m>2.0)";
    cs.addCut(new cutConst("(" + accstr1m_1 + "||" + accstr1m_2 +  "||" + accstr1m_3 + ")"));
    parvec.push_back(0);
    const std::string accstr2m_1 = "(TMath::Abs(Seta2m)<1.2&&Spt2m>3.5)";
    const std::string accstr2m_2 = "(TMath::Abs(Seta2m)>=1.2&&TMath::Abs(Seta2m)<1.6&&Sp2m>8-3.75*TMath::Abs(Seta2m))";
    const std::string accstr2m_3 = "(TMath::Abs(Seta2m)>=1.6&&TMath::Abs(Seta2m)<2.4&&Spt2m>2.0)";
    cs.addCut(new cutConst("(" + accstr2m_1 + "||" + accstr2m_2 +  "||" + accstr2m_3 + ")"));
    parvec.push_back(0);
    cs.addCut(new cutConst("(vrKs>1&&vrKs<35&&TMath::Abs(vzKs)<100)"));
    parvec.push_back(0);
}

void Cuts::selectCutAcceptance111215global50() // muons based on AN-11-417, see Journal #13 p66 for calculation of cuts
{
    const std::string accstr1m_1 = "(TMath::Abs(Seta1m)<1.2&&Spt1m>4.6-1.45*TMath::Abs(Seta1m))";
    const std::string accstr1m_2 = "(TMath::Abs(Seta1m)>=1.2&&TMath::Abs(Seta1m)<1.6&&Sp1m>8.59-4.17*TMath::Abs(Seta1m))";
    const std::string accstr1m_3 = "(TMath::Abs(Seta1m)>=1.6&&TMath::Abs(Seta1m)<2.4&&Spt1m>3.8-0.75*TMath::Abs(Seta1m))";
    cs.addCut(new cutConst("(" + accstr1m_1 + "||" + accstr1m_2 +  "||" + accstr1m_3 + ")"));
    parvec.push_back(0);
    const std::string accstr2m_1 = "(TMath::Abs(Seta2m)<1.2&&Spt2m>4.6-1.45*TMath::Abs(Seta2m))";
    const std::string accstr2m_2 = "(TMath::Abs(Seta2m)>=1.2&&TMath::Abs(Seta2m)<1.6&&Sp2m>8.59-4.17*TMath::Abs(Seta2m))";
    const std::string accstr2m_3 = "(TMath::Abs(Seta2m)>=1.6&&TMath::Abs(Seta2m)<2.4&&Spt2m>3.8-0.75*TMath::Abs(Seta2m))";
    cs.addCut(new cutConst("(" + accstr2m_1 + "||" + accstr2m_2 +  "||" + accstr2m_3 + ")"));
    parvec.push_back(0);
}

void Cuts::selectCutAcceptance111215_Ks() // cut on Ks decay vertex
{
    cs.addCut(new cutConst("(vrKs>1&&vrKs<35&&TMath::Abs(vzKs)<100)"));
    parvec.push_back(0);
}

void Cuts::selectCutAcceptance111215_L0() // cut on Ks decay vertex
{
    cs.addCut(new cutConst("(vrl0>1&&vrl0<35&&TMath::Abs(vzl0)<100)"));
    parvec.push_back(0);
}

// ask for barrel muons 
void Cuts::selectCutMuEta16()
{
    cs.addCut(new cutConst("TMath::Abs(Seta1m)<1.6&&TMath::Abs(Seta2m)<1.6"));
    parvec.push_back(0);
}

void Cuts::selectCutMuEta10()
{
    cs.addCut(new cutConst("TMath::Abs(Seta1m)<1.0&&TMath::Abs(Seta2m)<1.0"));
    parvec.push_back(0);
}

// trigger cuts
void Cuts::selectCutTrgL1_110120_0()
{
    cs.addCut(new cutConst("(L1TDMu0==1||L1TMu0==1)"));
    parvec.push_back(0);
}

void Cuts::selectCutTrgL1_110120_3()
{
    cs.addCut(new cutConst("(L1TDMu3==1||L1TMu3==1)"));
    parvec.push_back(0);
}

void Cuts::selectCutTrgHLT_110120()
{
    cs.addCut(new cutConst("(HLTqrk==1||HLTDMu3==1||HLTMu0jp==1||HLTMu0jpT==1||HLTMu3jp==1||HLTMu0jpT==1)"));
    parvec.push_back(0);
}

void Cuts::selectCutTrgHLT_110131()
{
    cs.addCut(new cutConst("(HLTMu0jp==1||HLTMu0jpT==1)"));
    parvec.push_back(0);
}

void Cuts::selectCutTrgMatch_110209()
{
    cs.addCut(new cutConst("HLTmatch==1&&HLTok==1"));
    parvec.push_back(0);
}

void Cuts::selectCutTrgMatch_2011_110413()
{
    cs.addCut(new cutConst("HLTqrk==1&&HLTmatch==1"));
    //cs.addCut(new cutConst("HLTmatch==1"));
    parvec.push_back(0);
}

void Cuts::selectCutTrgMatch_2011_110512()
{
    cs.addCut(new cutConst("HLTmatch==1"));
    //cs.addCut(new cutConst("HLTmatch==1"));
    parvec.push_back(0);
}

void Cuts::selectCutTrgMatch() { cs.addCut(new cutConst("HLTmatch==1")); parvec.push_back(0); }
void Cuts::selectCutTrgJpsi() { cs.addCut(new cutConst("HLTokJpsi==1")); parvec.push_back(0); }
void Cuts::selectCutTrgJpsiDispl() { cs.addCut(new cutConst("HLTokDisplJpsi==1")); parvec.push_back(0); }
void Cuts::selectCutTrgJpsiBarrel() { cs.addCut(new cutConst("HLTokBarrelJpsi==1")); parvec.push_back(0); }

// some special cuts
void Cuts::selectCutMlbWindow_110215()
{
    cs.addCut(new cutSymWindow("mlb",5.624));   parvec.push_back(0.022);
}

void Cuts::selectCutIso120229()
{
    cs.addCut(new cutBoundUpper("isoClostrk"));   parvec.push_back(21);
    cs.addCut(new cutBoundLower("isoDocatrk"));   parvec.push_back(0);
}

// cuts to select a certain bin in parameter space
void Cuts::selectCutPtLb_10_15()
{
    cs.addCut(new cutBoundLower("ptlb")); parvec.push_back(10.);
    cs.addCut(new cutBoundUpper("ptlb")); parvec.push_back(15.);
}
void Cuts::selectCutPtLb_15_20()
{
    cs.addCut(new cutBoundLower("ptlb")); parvec.push_back(15.);
    cs.addCut(new cutBoundUpper("ptlb")); parvec.push_back(20.);
}
void Cuts::selectCutPtLb_20_infty()
{
    cs.addCut(new cutBoundLower("ptlb")); parvec.push_back(20.);
}
void Cuts::selectCutYlb_00_05()
{
    cs.addCut(new cutConst("TMath::Abs(ylb)<0.5")); parvec.push_back(0);
}
void Cuts::selectCutYlb_05_12()
{
    cs.addCut(new cutConst("TMath::Abs(ylb)>=0.5")); parvec.push_back(0);
    cs.addCut(new cutConst("TMath::Abs(ylb)<1.2")); parvec.push_back(0);
}
void Cuts::selectCutYlb_12_22()
{
    cs.addCut(new cutConst("TMath::Abs(ylb)>=1.2")); parvec.push_back(0);
    cs.addCut(new cutConst("TMath::Abs(ylb)<2.2")); parvec.push_back(0);
}

