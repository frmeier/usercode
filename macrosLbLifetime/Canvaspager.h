#ifndef CANVASPAGER_H_GUARD
#define CANVASPAGER_H_GUARD
#include "utils.h"
#include "TCanvas.h"
#include <string>

class Canvaspager
{
    public:
	Canvaspager(TCanvas *c, std::string path, std::string filename, const int ncd)
	    : c_(c), path_(path), filename_(filename), ncd_(ncd), cdcounter_(0), pagecounter_(0) {};

	void cdNext()
	{
	    cdcounter_++;
	    if (cdcounter_ > ncd_)
	    {
		forceSave();
		pagecounter_++;
		cdcounter_=1;
		c_->Clear("D");
	    }
	    c_->cd(cdcounter_);
	};
	void cdNext(std::string name)
	{
	    cdNext();
	    name_ = name;
	};
	void cd() { c_->cd(cdcounter_); };
	void forceSave()
	{
	    const std::string filename = path_ + "/" + filename_ + toString(pagecounter_) + (name_.size() != 0 ? "_" + name_ : "");
	    c_->SaveAs((filename + ".pdf").c_str());
	    c_->SaveAs((filename + ".png").c_str());
	    lastSaved_ = filename;
	};
	int getNcd() { return ncd_; };
	std::string getLastSavedPng() { return lastSaved_ + ".png"; };
	std::string getCurPng() { return filename_ + toString(pagecounter_) + (name_.size() != 0 ? "_" + name_ : "") + ".png"; };
	std::string getLastSavedPdf() { return lastSaved_ + ".pdf"; };
	std::string getCurPdf() { return filename_ + toString(pagecounter_) + (name_.size() != 0 ? "_" + name_ : "") + ".pdf"; };

    private:
	TCanvas *c_;
	std::string path_, filename_;
	int ncd_;
	int cdcounter_;
	int pagecounter_;
	std::string lastSaved_;
	std::string name_;
};

#endif

