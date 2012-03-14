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

	void printCut()
	{
	};

	void selectCut(std::string name1, std::string name2 = "", std::string name3 = "", std::string name4 = "")
	{
	    init();
	    cout << "selectCut: ";
	    select1Cut(name1); cout << name1 << " ";
	    if (name2 != "") { select1Cut(name2); cout << name2 << " "; }
	    if (name3 != "") { select1Cut(name3); cout << name3 << " "; }
	    if (name4 != "") { select1Cut(name4); cout << name4 << " "; }
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
	void selectCut110525();
	void selectCut120309();
	void selectCut110203exp();
	void selectCut110203exp2();
	void selectCut110222exp();
	void selectCutB0_110628();
	void selectCutB0_110628exp();
	void selectCutAcceptance110120();
	void selectCutAcceptance110131();
	void selectCutAcceptance110210();
	void selectCutAcceptance110210B0();
	void selectCutAcceptance111215tracker50();
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

