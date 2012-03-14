#ifndef HTMLREPORT_GUARD_H
#define HTMLREPORT_GUARD_H

#include <string>
#include <fstream>

class HtmlReport
{
    public:
	HtmlReport() {};
	HtmlReport(std::string file, std::string title);
	~HtmlReport();

	void addH(std::string text, char level = '1', bool foldable = false);
	void addP(std::string text, bool foldable = false);

	void beginTable();
	void endTable();
	void addTableImage(std::string file, std::string caption, std::string href = "");
	void addTableCell(std::string text);
	void beginTableRow();
	void endTableRow();
	void addTableRow(std::string t1, std::string t2 = "", std::string t3 = "", std::string t4 = "");

	void beginDiv();
	void endDiv();
	void flushDiv();
	void addImage(std::string file, std::string caption);
	void addImage(std::string file, std::string caption, std::string href);

    private:
	std::string file_;
	std::ofstream out_;

	bool isTable_;
	bool isTableRow_;
	unsigned int divCounter_;
	bool isCollapsable_;

	void initPage(std::string title);
	void finalizePage();
};

#endif

