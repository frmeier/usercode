#ifndef JSONRUNLIST_H_GUARD
#define JSONRUNLIST_H_GUARD

#include <set>
#include <map>
#include <iostream>
#include <ostream>

class JsonRunList
{
    public:
	typedef std::map<int, std::set<int> > runlistType;
	JsonRunList() {};
	void insert(int run, int ls)
	{
	    runlist_[run].insert(ls);
	};
	int getSize()
	{
	    int ret(0);
	    for (runlistType::const_iterator itmap = runlist_.begin(); itmap!=runlist_.end(); itmap++)
	    {
		ret += itmap->second.size();
	    }
	    return ret;
	};
	friend std::ostream& operator<<(std::ostream& os, const JsonRunList &jsr)
	{
	    bool moreRun(false);
	    os << "{";
	    for (runlistType::const_iterator itmap = jsr.runlist_.begin(); itmap!=jsr.runlist_.end(); itmap++)
	    {
		if (moreRun) os << "," << std::endl;
		os << "\"" << itmap->first << "\" : [";
		bool moreLs(false);
		bool loSet(false);
		int lastVal(0);
		for (std::set<int>::const_iterator itset = itmap->second.begin(); itset != itmap->second.end(); itset++)
		{
		    if (!loSet)
		    {
			if (moreLs) os << " , ";
			lastVal = (*itset);
			os << "[" << (*itset) << ",";
			loSet = true;
		    }
		    else
		    {
			if (lastVal+1 != (*itset)) 
			{
			    os << lastVal << "]";
			    os << ",[" << (*itset) << ",";
			}
			lastVal = (*itset);
		    }
		    moreLs =true;
		}
		os << lastVal << "]]";
		moreRun = true;
	    }
	    os << "}" << std::endl;
	    return os;
	};

    //private:
	runlistType runlist_;
};

#endif

