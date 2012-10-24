#ifndef CUTS_H_GUARD
#define CUTS_H_GUARD

#include "Cut.h"

class Cuts
{
    public:
	typedef cutBase::valueType valueType;
	typedef std::vector<cutBase::valueType> parVecType;
	cutSet cs;
	parVecType parvec;
	cutBase::cutnameType getCut() const { return cs.getCut(parvec); };
	cutBase::cutnameType getOneCut(parVecType::size_type i)
	{
	    if (i<cs.size()) // i>=0 because of unsignedness
		return cs.getOneCut(i, parvec[i]);
	    else
		throw ("Out of bounds");
	};
	cutBase::cutnameType getOneCut(cutBase::cutnameType name)
	{
	    return getOneCut(cs.getCutPos(name));
	};

	cutBase::valueType getOneCutValue(cutBase::cutnameType name)
	{
	    const parVecType::size_type i = cs.getCutPos(name);
	    if (i<parvec.size()) // i>=0 because of unsignedness
		return parvec[i];
	    else
		throw ("Out of bounds");
	};

	void removeOneCut(cutBase::cutnameType name)
	{
	    const parVecType::size_type i = cs.getCutPos(name);
	    if (i<parvec.size()) // i>=0 because of unsignedness
	    {
		cs.erase(cs.begin() + i);
		parvec.erase(parvec.begin() + i);
	    }
	    else
		throw ("Out of bounds");
	};

	void printCutLaTeX_values(std::string prefix = "") { cs.printCutLaTeX_values(parvec, prefix); };
	
	cutBase::cutnameType getOneCut(cutBase::cutnameType name, cutBase::valueType v)
	{
	    return cs.getOneCut(cs.getCutPos(name),v);
	};
	cutBase::cutnameType getCutExceptOne(cutBase::cutnameType name)
	{
	    return cs.getCutExceptOne(parvec, cs.getCutPos(name));
	};
	cutBase::cutnameType getCutChangeOne(cutBase::cutnameType name, cutBase::valueType v)
	{
	    return getCutExceptOne(name) + "&&" + getOneCut(name, v);
	};

	void printCut()
	{
	};

	void selectCut(std::string name1, std::string name2 = "", std::string name3 = "", std::string name4 = "",
	       	std::string name5 = "", std::string name6 = "", std::string name7 = "")
	{
	    init();
	    cout << "selectCut: ";
	    select1Cut(name1); cout << name1 << " ";
	    if (name2 != "") { select1Cut(name2); cout << name2 << " "; }
	    if (name3 != "") { select1Cut(name3); cout << name3 << " "; }
	    if (name4 != "") { select1Cut(name4); cout << name4 << " "; }
	    if (name5 != "") { select1Cut(name5); cout << name5 << " "; }
	    if (name6 != "") { select1Cut(name6); cout << name6 << " "; }
	    if (name7 != "") { select1Cut(name7); cout << name7 << " "; }
	    cout << endl;
	};

	void init()
	{
	    cs.clear();
	    parvec.clear();
	};

	Cuts operator+=(const Cuts &c)
	{
	    cs+=c.cs;
	    parvec.insert(parvec.end(), c.parvec.begin(), c.parvec.end());
	    return *this;
	};

	Cuts operator+(const Cuts &c)
	{
	    Cuts ret = *this;
	    ret+=c;
	    return ret;
	};

    private:
	void select1Cut(std::string name);
	void selectCutTrivial();
	void selectCut101208();
	void selectCut110131();
	void selectCut110203();
	void selectCut110214();
	void selectCut110303();
	void selectCut110303exp();
	void selectCut110413();
	void selectCutLb110525();
	void selectCutLb120309();
	void selectCutLb120315();
	void selectCutLb120417();
	void selectCutLb120420();
	void selectCutLb120523();
	void selectCutLb120711();
	void selectCutLb120713();
	void selectCutLb120905();
	void selectCutLb120610();
	void selectCut110203exp();
	void selectCut110203exp2();
	void selectCut110222exp();
	void selectCutB0_110628();
	void selectCutB0_110628exp();
	void selectCutB0_120315();
	void selectCutB0_120414();
	void selectCutB0_120416();
	void selectCutB0_120420();
	void selectCutB0_120523();
	void selectCutB0_120711();
	void selectCutB0_120905();
	void selectCutAcceptance110120();
	void selectCutAcceptance110131();
	void selectCutAcceptance110210();
	void selectCutAcceptance110210B0();
	void selectCutAcceptance111215tracker50();
	void selectCutAcceptance120523tracker50();
	void selectCutAcceptance111215global50();
	void selectCutAcceptance111215_Ks();
	void selectCutAcceptance111215_L0();
	void selectCutMuEta16();
	void selectCutMuEta10();
	void selectCutTrgL1_110120_0();
	void selectCutTrgL1_110120_3();
	void selectCutTrgHLT_110120();
	void selectCutTrgHLT_110131();
	void selectCutTrgMatch_110209();
	void selectCutTrgMatch_2011_110413();
	void selectCutTrgMatch_2011_110512();
	void selectCutTrgMatch();
	void selectCutTrgJpsi();
	void selectCutTrgJpsiDispl();
	void selectCutTrgJpsiBarrel();
	void selectCutTrgJpsiDisplMCPseudo();
	void selectCutTrgJpsiBarrelMCPseudo();
	void selectCutMlbWindow_110215();
	void selectCutIso120229();
	void selectCutPtLb_10_15();
	void selectCutPtLb_15_20();
	void selectCutPtLb_20_infty();
	void selectCutYlb_00_05();
	void selectCutYlb_05_12();
	void selectCutYlb_12_22();
};

#endif

