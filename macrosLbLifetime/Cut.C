#include <string>
#include <sstream>
#include <vector>
#include <utility>
#include <iostream>
#include <stdexcept>
#include "Math/SVector.h"
#include "Cut.h"
#include "utils.h"

cutBase::cutnameType cutSet::getOneCut(sizeType i, cutBase::valueType v) const
{
    if (i<cutvec_.size()) // i>=0 because of unsignedness
	return cutvec_[i]->getCut(v);
    else
	throw ("Out of bounds");
};

void cutTest()
{
    cutSet cs;
    cs.addCut(new cutGT2Var("ptpr","ptpi"));
    cs.addCut(new cutSymWindow("mjp",3));
    cs.addCut(new cutBoundLower("ptpr"));

    std::vector<cutBase::valueType> parvec;
    parvec.push_back(0);
    parvec.push_back(.1);
    parvec.push_back(8);
    cout << cs.getCut(parvec) << endl;
    cout << cs.getCutInvertN(parvec,1) << endl;
}

