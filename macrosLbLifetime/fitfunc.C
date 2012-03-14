#include "TF1.h"
#include "TH1F.h"

class fitfuncBase
{
    public:
	typedef double valueType;
	virtual void doFit(TH1F *h) = 0;
	virtual void reset() = 0;
	bool isValid() const { return isValid_; };
	valueType getSig() const { return sig_; };
	valueType getBgr() const { return bgr_; };
	valueType getSoverSqrtSB() const { return sig_/sqrt(sig_+bgr_); };
	valueType getChi2() const { return chi2_; };
	virtual valueType getParameter(const int i) const = 0;

    protected:
	bool isValid_;
	valueType chi2_;
	valueType min_, max_;
	valueType sig_, bgr_;
};

class fitfuncGausBgrFlat : protected fitfuncBase
{
    public:
	fitfuncGausBgrFlat() : nSigmas_(2.0)
	{
	    init();
	};
	fitfuncGausBgrFlat( valueType min, valueType max, 
		valueType pConst, valueType pSlope, valueType pScale, valueType pMean, valueType pSigma,
		valueType nSigmas) :
	    pConst_(pConst), pSlope_(pSlope), pScale_(pScale), pMean_(pMean), pSigma_(pSigma), nSigmas_(nSigmas)
	{
	    min_ = min;
	    max_ = max;
	    init();
	};
	//~fitfuncGausBgrFlat()
	//{
	    //f1_->delete;
	    //f2_->delete;
	//}

	void doFit(TH1F *h)
	{
	    // do the fit
	    h->Fit(f1_,"R");
	    // extract results
	    const valueType binwidth = h->GetBinWidth(1); // assuming equally distributed bins, may break!
	    const valueType scaleResult = f1_->GetParameter(2);
	    const valueType meanResult  = f1_->GetParameter(3);
	    const valueType sigmaResult = f1_->GetParameter(4);
	    const valueType lowerBound = meanResult-nSigmas_*sigmaResult;
	    const valueType upperBound = meanResult+nSigmas_*sigmaResult;
	    const valueType SB = f1_->Integral(lowerBound,upperBound)/binwidth;
	    // function for the signal
	    f2_->SetParameters(scaleResult,meanResult,sigmaResult);
	    const valueType S = f2_->Integral(lowerBound,upperBound)/binwidth;
	    // set values
	    sig_ = S;
	    bgr_ = SB-S;
	    chi2_ = f1_->GetChisquare();
	    isValid_ = true;
	    cout << "lowerBound: " << lowerBound << " upperBound: " << upperBound;
	    cout << " S: " << S << " SB: " << SB << endl;
	};
	void init()
	{
	    //if(f1_ != 0) delete f1_; 
	    //if(f2_ != 0) delete f2_;
	    // fit function: Gaussian plus flat bgr
	    f1_ = new TF1("f1", "[0] + [1]*x + gaus(2)", min_, max_);
	    f1_->SetParNames("const","slope","scale","mean","sigma");
	    // function for signal
	    f2_ = new TF1("f2", "gaus", min_, max_);
	    f2_->SetParNames("scale","mean","sigma");
	};
	void setInitValues(valueType min, valueType max, 
		valueType pConst, valueType pSlope, valueType pScale, valueType pMean, valueType pSigma,
		valueType nSigmas)
	{
	    min_ = min;
	    max_ = max;
	    pConst_ = pConst;
	    pSlope_ = pSlope;
	    pScale_ = pScale;
	    pMean_  = pMean;
	    pSigma_ = pSigma;
	    nSigmas_ = nSigmas;
	    init();
	};

	void reset()
	{
	    f1_->SetParameters(pConst_, pSlope_, pScale_, pMean_, pSigma_);
	};

	valueType getParameter(const int i) const { return f1_->GetParameter(i); } ;
	valueType get_pConst() { return f1_->GetParameter(0); };
	valueType get_pSlope() { return f1_->GetParameter(1); };
	valueType get_pScale() { return f1_->GetParameter(2); };
	valueType get_pMean()  { return f1_->GetParameter(3); };
	valueType get_pSigma() { return f1_->GetParameter(4); };

    private:
	// parameters of the fit function
	valueType pConst_;
	valueType pSlope_;
	valueType pScale_;
	valueType pMean_;
	valueType pSigma_;
	// width of gaussian used for S and B estimation
	valueType nSigmas_;
	// functions
	TF1 *f1_, *f2_;
};

class fitfuncGausBgrExp : protected fitfuncBase
{
    public:
	fitfuncGausBgrExp() : nSigmas_(2.0)
	{
	    init();
	};
	fitfuncGausBgrExp( valueType min, valueType max, 
		valueType pConst, valueType pDecay, valueType pScale, valueType pMean, valueType pSigma,
		valueType nSigmas) :
	    pConst_(pConst), pDecay_(pDecay), pScale_(pScale), pMean_(pMean), pSigma_(pSigma), nSigmas_(nSigmas)
	{
	    min_ = min;
	    max_ = max;
	    init();
	};
	~fitfuncGausBgrExp()
	{
	    delete f1_;
	    delete f2_;
	}

	void doFit(TH1F *h)
	{
	    // do the fit
	    h->Fit(f1_,"R");
	    // extract results
	    const valueType binwidth = h->GetBinWidth(1); // assuming equally distributed bins, may break!
	    const valueType scaleResult = f1_->GetParameter(2);
	    const valueType meanResult  = f1_->GetParameter(3);
	    const valueType sigmaResult = f1_->GetParameter(4);
	    const valueType lowerBound = meanResult-nSigmas_*sigmaResult;
	    const valueType upperBound = meanResult+nSigmas_*sigmaResult;
	    const valueType SB = f1_->Integral(lowerBound,upperBound)/binwidth;
	    // function for the signal
	    f2_->SetParameters(scaleResult,meanResult,sigmaResult);
	    const valueType S = f2_->Integral(lowerBound,upperBound)/binwidth;
	    // set values
	    sig_ = S;
	    bgr_ = SB-S;
	    chi2_ = f1_->GetChisquare();
	    isValid_ = true;
	    cout << "lowerBound: " << lowerBound << " upperBound: " << upperBound;
	    cout << " S: " << S << " SB: " << SB << endl;
	};
	void init()
	{
	    //if(f1_ != 0) delete f1_;
	    //if(f2_ != 0) delete f2_;
	    // fit function: Gaussian plus exponential bgr
	    f1_ = new TF1("f1", "[0] * exp([1]*x) + gaus(2)", min_, max_);
	    f1_->SetParNames("const","decay","scale","mean","sigma");
	    // function for signal
	    f2_ = new TF1("f2", "gaus", min_, max_);
	    f2_->SetParNames("scale","mean","sigma");
	};
	void setInitValues(valueType min, valueType max, 
		valueType pConst, valueType pDecay, valueType pScale, valueType pMean, valueType pSigma,
		valueType nSigmas)
	{
	    min_ = min;
	    max_ = max;
	    pConst_ = pConst;
	    pDecay_ = pDecay;
	    pScale_ = pScale;
	    pMean_  = pMean;
	    pSigma_ = pSigma;
	    nSigmas_ = nSigmas;
	    init();
	};

	void reset()
	{
	    f1_->SetParameters(pConst_, pDecay_, pScale_, pMean_, pSigma_);
	};

	valueType getParameter(const int i) const { return f1_->GetParameter(i); } ;
	valueType get_pConst() { return f1_->GetParameter(0); };
	valueType get_pDecay() { return f1_->GetParameter(1); };
	valueType get_pScale() { return f1_->GetParameter(2); };
	valueType get_pMean()  { return f1_->GetParameter(3); };
	valueType get_pSigma() { return f1_->GetParameter(4); };

    private:
	// parameters of the fit function
	valueType pConst_;
	valueType pDecay_;
	valueType pScale_;
	valueType pMean_;
	valueType pSigma_;
	// width of gaussian used for S and B estimation
	valueType nSigmas_;
	// functions
	TF1 *f1_, *f2_;
};

class fitfuncGausBgrExpConst : protected fitfuncBase
{
    public:
	fitfuncGausBgrExpConst() : nSigmas_(2.0)
	{
	    init();
	};
	fitfuncGausBgrExpConst( valueType min, valueType max, 
		valueType pExpscale, valueType pDecay, valueType pScale, valueType pMean, valueType pSigma, valueType pConst, 
		valueType nSigmas) :
	    pExpscale_(pExpscale), pDecay_(pDecay), pScale_(pScale), pMean_(pMean), pSigma_(pSigma), pConst_(pConst),
	    nSigmas_(nSigmas)
	{
	    min_ = min;
	    max_ = max;
	    init();
	};
	~fitfuncGausBgrExpConst()
	{
	    delete f1_;
	    delete f2_;
	}

	void doFit(TH1F *h)
	{
	    // do the fit
	    h->Fit(f1_,"R");
	    // extract results
	    const valueType binwidth = h->GetBinWidth(1); // assuming equally distributed bins, may break!
	    const valueType scaleResult = f1_->GetParameter(2);
	    const valueType meanResult  = f1_->GetParameter(3);
	    const valueType sigmaResult = f1_->GetParameter(4);
	    const valueType lowerBound = meanResult-nSigmas_*sigmaResult;
	    const valueType upperBound = meanResult+nSigmas_*sigmaResult;
	    const valueType SB = f1_->Integral(lowerBound,upperBound)/binwidth;
	    // function for the signal
	    f2_->SetParameters(scaleResult,meanResult,sigmaResult);
	    const valueType S = f2_->Integral(lowerBound,upperBound)/binwidth;
	    // set values
	    sig_ = S;
	    bgr_ = SB-S;
	    chi2_ = f1_->GetChisquare();
	    isValid_ = true;
	    cout << "lowerBound: " << lowerBound << " upperBound: " << upperBound;
	    cout << " S: " << S << " SB: " << SB << endl;
	};
	void init()
	{
	    if(f1_ != 0) delete f1_;
	    if(f2_ != 0) delete f2_;
	    // fit function: Gaussian plus exponential bgr
	    f1_ = new TF1("f1", "[0] * exp([1]*x) + gaus(2) + [5]", min_, max_);
	    f1_->SetParNames("expscale","decay","scale","mean","sigma","const");
	    // function for signal
	    f2_ = new TF1("f2", "gaus", min_, max_);
	    f2_->SetParNames("scale","mean","sigma");
	};
	void setInitValues(valueType min, valueType max, 
		valueType pExpscale, valueType pDecay, valueType pScale, valueType pMean, valueType pSigma, valueType pConst,
		valueType nSigmas)
	{
	    min_ = min;
	    max_ = max;
	    pExpscale_ = pExpscale;
	    pDecay_ = pDecay;
	    pScale_ = pScale;
	    pMean_  = pMean;
	    pSigma_ = pSigma;
	    pConst_ = pConst;
	    nSigmas_ = nSigmas;
	    init();
	};

	void reset()
	{
	    f1_->SetParameters(pExpscale_, pDecay_, pScale_, pMean_, pSigma_, pConst_);
	    // set some limits
	    //f1_->SetParLimits(3,0,5*max_); // mean must be positive
	    f1_->SetParLimits(4,0,(max_-min_)*100); // sigma must be positive and not wider than the plot
	    //f1_->SetParLimits(5,0,100); // const must be positive
	};

	valueType getParameter(const int i) const { return f1_->GetParameter(i); } ;
	valueType get_pExpscale() { return f1_->GetParameter(0); };
	valueType get_pDecay() { return f1_->GetParameter(1); };
	valueType get_pScale() { return f1_->GetParameter(2); };
	valueType get_pMean()  { return f1_->GetParameter(3); };
	valueType get_pSigma() { return f1_->GetParameter(4); };
	valueType get_pConst() { return f1_->GetParameter(5); };

    private:
	// parameters of the fit function
	valueType pExpscale_;
	valueType pDecay_;
	valueType pScale_;
	valueType pMean_;
	valueType pSigma_;
	valueType pConst_;
	// width of gaussian used for S and B estimation
	valueType nSigmas_;
	// functions
	TF1 *f1_, *f2_;
};

