#ifndef CUT_H_GUARD
#define CUT_H_GUARD

#include <string>
#include <sstream>
#include <vector>
#include <utility>
#include <iostream>
#include <stdexcept>
#include "Math/SVector.h"
#include "utils.h"

class cutBase
{
    public:
	typedef double valueType;
	typedef std::string cutnameType;

	cutBase();
	cutBase(cutnameType s) : cutname_(s) {};
	virtual cutnameType getCut(valueType v) const = 0;
	virtual valueType getCutValue(valueType v) const = 0;
	virtual cutnameType getCutInverted(valueType v) const = 0;
	cutnameType getCutName() const { return cutname_; } ;
	virtual cutnameType getCutTitle() const = 0 ;
	virtual cutnameType getCutClassName() const = 0 ;
	virtual cutnameType getCutLaTeX(valueType v) const = 0;

    protected:
	cutnameType cutname_;
};

class cutConst : public cutBase
{
    public:
	cutConst(cutnameType s) : cutBase(s) {};
	cutnameType getCut(valueType v = 0) const { (void)v; return cutname_; };
	valueType getCutValue(valueType v = 0) const { (void)v; return 0; };
	cutnameType getCutInverted(valueType v = 0) const { (void)v; return cutname_; };
	cutnameType getCutTitle() const { return cutname_; };
	cutnameType getCutClassName() const { return "cutConst"; };
	cutnameType getCutLaTeX(valueType v = 0) const { (void)v; return cutname_; };
};

class cutGT2Var : public cutBase
{
    public:
	cutGT2Var(cutnameType lhs, cutnameType rhs) : cutBase(lhs), rhs_(rhs) {};
	cutnameType getCut(valueType v = 0) const { (void)v; return cutname_ + ">" + rhs_; };
	valueType getCutValue(valueType v = 0) const { (void)v; return 0; };
	cutnameType getCutInverted(valueType v = 0) const { (void)v; return cutname_ + "<=" + rhs_; };
	cutnameType getCutTitle() const { return cutname_ + ">" + rhs_; };
	cutnameType getCutClassName() const { return "cutGT2Var"; };
	cutnameType getCutLaTeX(valueType v = 0) const { (void)v; return cutname_ + ">" + rhs_; };
    private:
	cutnameType rhs_;
};

class cutLT2Var : public cutBase
{
    public:
	cutLT2Var(cutnameType lhs, cutnameType rhs) : cutBase(lhs), rhs_(rhs) {};
	cutnameType getCut(valueType v = 0) const { (void)v; return cutname_ + "<" + rhs_; };
	valueType getCutValue(valueType v = 0) const { (void)v; return 0; };
	cutnameType getCutInverted(valueType v = 0) const { (void)v; return cutname_ + ">=" + rhs_; };
	cutnameType getCutTitle() const { return cutname_ + "<" + rhs_; };
	cutnameType getCutClassName() const { return "cutLT2Var"; };
	cutnameType getCutLaTeX(valueType v = 0) const { (void)v; return cutname_ + "<" + rhs_; };
    private:
	cutnameType rhs_;
};

class cutEqual : public cutBase
{
    public:
	cutEqual(cutnameType s) : cutBase(s) {};
	cutnameType getCut(valueType v) const { return cutname_ + "==" + toString(v); };
	valueType getCutValue(valueType v = 0) const { return v; };
	cutnameType getCutInverted(valueType v) const { return cutname_ + "!=" + toString(v); };
	cutnameType getCutTitle() const { return cutname_ + "=="; };
	cutnameType getCutClassName() const { return "cutEqual"; };
	cutnameType getCutLaTeX(valueType v = 0) const { return "= " + toString(v); };
};

class cutBoundLower : public cutBase
{
    public:
	cutBoundLower(cutnameType s) : cutBase(s) {};
	cutnameType getCut(valueType v) const { return cutname_ + ">" + toString(v); }
	valueType getCutValue(valueType v = 0) const { return v; };
	cutnameType getCutInverted(valueType v) const { return cutname_ + "<=" + toString(v); }
	cutnameType getCutTitle() const { return cutname_ + ">"; };
	cutnameType getCutClassName() const { return "cutBoundLower"; };
	cutnameType getCutLaTeX(valueType v = 0) const { return "> " + toString(v); };
};

class cutBoundUpper : public cutBase
{
    public:
	cutBoundUpper(cutnameType s) : cutBase(s) {};
	cutnameType getCut(valueType v) const { return cutname_ + "<" + toString(v); }
	valueType getCutValue(valueType v = 0) const { return v; };
	cutnameType getCutInverted(valueType v) const { return cutname_ + ">=" + toString(v); }
	cutnameType getCutTitle() const { return cutname_ + "<"; };
	cutnameType getCutClassName() const { return "cutBoundUpper"; };
	cutnameType getCutLaTeX(valueType v = 0) const { return "< " + toString(v); };
};

class cutSymWindow : public cutBase
{
    public:
	cutSymWindow(cutnameType s, valueType v) : cutBase(s), center_(v) {};
	cutnameType getCut(valueType v) const
	{
	    return cutname_ + ">" + toString(center_-v) + "&&" + cutname_ + "<" + toString(center_+v);
	}
	valueType getCutValue(valueType v = 0) const { return center_ + v; };
	cutnameType getCutInverted(valueType v) const
	{
	    return "(" + cutname_ + "<=" + toString(center_-v) + "||" + cutname_ + ">=" + toString(center_+v) + ")";
	}
	cutnameType getCutTitle() const { return cutname_ + " window around " + toString(center_); };
	cutnameType getCutClassName() const { return "cutSymWindow"; };
	cutnameType getCutLaTeX(valueType v = 0) const { return "\\pm " + toString(v); };
    private:
	valueType center_;
};

class cutSymWindowVeto : public cutBase
{
    public:
	cutSymWindowVeto(cutnameType s, valueType v) : cutBase(s), center_(v) {};
	cutnameType getCut(valueType v) const
	{
	    return "(" + cutname_ + "<" + toString(center_-v) + "||" + cutname_ + ">" + toString(center_+v) + ")";
	}
	valueType getCutValue(valueType v = 0) const { return center_ + v; };
	cutnameType getCutInverted(valueType v) const
	{
	    return cutname_ + ">=" + toString(center_-v) + "&&" + cutname_ + "<=" + toString(center_+v);
	}
	cutnameType getCutTitle() const { return cutname_ + " veto window around " + toString(center_); };
	cutnameType getCutClassName() const { return "cutSymWindowVeto"; };
	cutnameType getCutLaTeX(valueType v = 0) const { return "\\pm " + toString(v); };
    private:
	valueType center_;
};

class cutRatioWindow : public cutBase
{   // ratio window: center = 1 and value = 2 gives a window from 1/2 to 2; value 3 from 1/3 to 3 etc.
    // center = 2 and value = 2: 1..4
    public:
	cutRatioWindow(cutnameType s, valueType v) : cutBase(s), center_(v) {};
	cutnameType getCut(valueType v) const
	{
	    return cutname_ + ">" + toString(center_/v) + "&&" + cutname_ + "<" + toString(center_*v);
	}
	valueType getCutValue(valueType v = 0) const { return center_ + v; };
	cutnameType getCutInverted(valueType v) const
	{
	    return "(" + cutname_ + "<=" + toString(center_/v) + "||" + cutname_ + ">=" + toString(center_*v) + ")";
	}
	cutnameType getCutTitle() const { return cutname_ + " ratio window around " + toString(center_); };
	cutnameType getCutClassName() const { return "cutRatioWindow"; };
	cutnameType getCutLaTeX(valueType v = 0) const { return "\\frac{"+toString(center_)+"}{"+ toString(v)+"}\\dots "+toString(center_)+"\\cdot"+toString(v); };
    private:
	valueType center_;
};

class cutBitcheck : public cutBase
{
    public:
	cutBitcheck(cutnameType s): cutBase(s) {};
	cutnameType getCut(valueType v) const
	{
	    return "(" + cutname_ + "&" + toString(v) + ")==" + toString(v);
	};
	valueType getCutValue(valueType v = 0) const { (void)v; return 0; };
	cutnameType getCutInverted(valueType v) const
	{
	    return "(" + cutname_ + "&" + toString(v) + ")!=" + toString(v);
	};
	cutnameType getCutTitle() const { return "bitcheck for " + cutname_; };
	cutnameType getCutClassName() const { return "cutBitcheck"; };
	cutnameType getCutLaTeX(valueType v = 0) const { return "bit(s) "+toString(v)+" set"; };
};

class cutSet
{
    public:
	typedef std::vector<cutBase *> cutvecType;
	typedef cutvecType::iterator iterator;
	typedef cutvecType::size_type sizeType;
	void addCut(cutBase *cut) { cutvec_.push_back(cut); };
	void clear() { cutvec_.clear(); };
	sizeType size() const { return cutvec_.size(); };
	
	cutSet operator+=(const cutSet cs)
	{
	    cutvec_.insert(cutvec_.end(), cs.cutvec_.begin(), cs.cutvec_.end());
	    return *this;
	};

	cutBase::cutnameType getOneCut(sizeType i, cutBase::valueType v) const;
	cutBase::valueType getOneCutValue(sizeType i, cutBase::valueType v) const
	{
	    if (i<cutvec_.size()) // i>=0 because of unsignedness
		return cutvec_[i]->getCutValue(v);
	    else
		throw ("Out of bounds");
	}
	cutBase::cutnameType getCutName(sizeType i) const { return cutvec_[i]->getCutName(); };
	cutBase::cutnameType getCutClassName(sizeType i) const { return cutvec_[i]->getCutClassName(); };
	sizeType getCutPos(cutBase::cutnameType cutname) const
	{
	    sizeType ret(0);
	    bool found(0);
	    for (sizeType i = 0; i!= cutvec_.size(); i++)
	    {
		if (cutvec_[i]->getCutName() == cutname)
		{
		    ret = i;
		    found = true;
		}
	    }
	    if (!found)
	    {
		cout << "Cutname " << cutname << " not found. Aborting" << endl;;
		throw std::invalid_argument(("Cutname " + cutname + " not found").c_str());
	    }
	    return ret;
       	};
	cutBase::cutnameType getCutTitle(sizeType i) const { return cutvec_[i]->getCutTitle(); };
	iterator begin() { return cutvec_.begin(); };
	iterator end() { return cutvec_.end(); };
	iterator erase(iterator it) { return cutvec_.erase(it); };
	// get the whole cut string
	cutBase::cutnameType getCut(std::vector<cutBase::valueType> pars) const
	{
	    if (pars.size() > 0)
		return getCut(&pars[0], pars.size());
	    else
		return "";
	};
	cutBase::cutnameType getCut(cutBase::valueType pars[], sizeType n) const
	{
	    if(n>cutvec_.size()) throw std::length_error("Too many parameters given");
	    cutBase::cutnameType ret("");
	    ret = cutvec_[0]->getCut(pars[0]);
	    for(sizeType i=1; i<n; i++)
	    {
		ret += "&&" + cutvec_[i]->getCut(pars[i]);
	    }
	    return ret;
	};
	// print the whole cut to cout as LaTeX code suitable for a values file
	void printCutLaTeX_values(std::vector<cutBase::valueType> pars, std::string prefix = "") const
	{
	    const unsigned int n = pars.size();
	    if(n>cutvec_.size()) throw std::length_error("Too many parameters given");
	    for(sizeType i=0; i<n; i++)
	    {
		cout << "\\vdef{" << prefix << cutvec_[i]->getCutName() << "}{\\ensuremath{" << cutvec_[i]->getCutLaTeX(pars[i]) << "}}" << endl;
		// \vdef{cuts:anal:probpr}{\ensuremath{>0.02}}
	    }
	};
	// get the whole cutstring with one cut explicitely inverted
	cutBase::cutnameType getCutInvertN(std::vector<cutBase::valueType> pars, sizeType invN) const
	{
	    return getCutInvertN(&pars[0], pars.size(), invN);
	};
	cutBase::cutnameType getCutInvertN(cutBase::valueType pars[], sizeType n, sizeType invN) const
	{
	    if(n>cutvec_.size()) throw std::length_error("Too many parameters given");
	    cutBase::cutnameType ret("");
	    if (invN == 0)
		ret = cutvec_[0]->getCutInverted(pars[0]);
	    else
		ret = cutvec_[0]->getCut(pars[0]);
	    for(sizeType i=1; i<n; i++)
	    {
		if (invN == i)
		    ret += "&&" + cutvec_[i]->getCutInverted(pars[i]);
		else
		    ret += "&&" + cutvec_[i]->getCut(pars[i]);
	    }
	    return ret;
	};
	// get the whole cutstring except one
	cutBase::cutnameType getCutExceptOne(std::vector<cutBase::valueType> pars, sizeType excI) const
	{
	    return getCutExceptOne(&pars[0], pars.size(), excI);
	};
	cutBase::cutnameType getCutExceptOne(cutBase::valueType pars[], sizeType n, sizeType excI) const
	{
	    if(n>cutvec_.size()) throw std::length_error("Too many parameters given");
	    cutBase::cutnameType ret("");
	    cutBase::cutnameType spacer("");
	    if (excI != 0)
	    {
		ret = cutvec_[0]->getCut(pars[0]);
		spacer = "&&";
	    }
	    for(sizeType i=1; i<n; i++)
	    {
		if (excI != i)
		{
		    ret += spacer + cutvec_[i]->getCut(pars[i]);
		    spacer = "&&";
		}
	    }
	    return ret;
	};
    private:
	cutvecType cutvec_;
};

#endif

