#include <string>
#include <fstream>
#include "HtmlReport.h"

HtmlReport::HtmlReport(std::string file, std::string title) : isTable_(false), isTableRow_(false), divCounter_(0), isCollapsable_(false)
{
    file_ = file;
    out_.open(file.c_str());
    initPage(title);
}

HtmlReport::~HtmlReport()
{
    finalizePage();
    out_.close();
}

void HtmlReport::initPage(std::string title)
{
    out_ << "<html>" << endl;
    out_ << "<head>";
    out_ << "<title>" << title << "</title>" << endl;
    out_ << "<style type=\"text/css\">" << endl;
    out_ << "@import \"../domcollapse.css\";" << endl;
    out_ << "</style>" << endl;
    out_ << "<script type=\"text/javascript\" src=\"../domcollapse.js\"></script>" << endl;
    out_ << "</head>" << endl;
    out_ << "<link rel=\"stylesheet\" type=\"text/css\" href=\"../frank.css\">";
    out_ << "<body>" << endl;
    return;
}

void HtmlReport::finalizePage()
{
    if (isTable_) endTable();
    flushDiv();
    out_ << "</body>" << endl;
    out_ << "</html>" << endl;
    return;
}

// =======================================================================================
void HtmlReport::beginTable()
{
    if (isTable_) endTable();
    beginDiv();
    out_ << "<table>" << endl;
    isTable_ = true;
    return;
}

void HtmlReport::endTable()
{
    if (!isTable_) return;
    if (isTableRow_)
    {
	out_ << "</tr>" << endl;
	isTable_ = false;
    }
    out_ << "</table>" << endl;
    isTable_ = false;
    endDiv();
    return;
}

// ---------------------------------------------------------------------------------------
void HtmlReport::beginTableRow()
{
    if (!isTable_) beginTable();
    if (isTableRow_) endTableRow();
    out_ << "<tr>" << endl;
    isTableRow_ = true;
    return;
}

void HtmlReport::endTableRow()
{
    if (!isTable_) return;
    if (isTableRow_) out_ << "</tr>" << endl;
    endDiv();
    isTableRow_ = false;
}

void HtmlReport::addTableRow(std::string t1, std::string t2, std::string t3, std::string t4)
{
    if (isTableRow_) endTableRow();
    addTableCell(t1);
    if (t2.size() != 0 || t3.size() != 0 || t4.size() !=0) addTableCell(t2);
    if (t3.size() != 0 || t4.size() !=0) addTableCell(t3);
    if (t4.size() !=0) addTableCell(t4);
    endTableRow();
}

// ---------------------------------------------------------------------------------------
void HtmlReport::addTableCell(std::string text)
{
    if (!isTable_) beginTable();
    if (!isTableRow_) beginTableRow();
    out_ << "<td>" << text << "</td>";
}

// ---------------------------------------------------------------------------------------
void HtmlReport::addTableImage(std::string file, std::string caption, std::string href)
{
    if (!isTable_) beginTable();
    if (!isTableRow_) beginTableRow();
    out_ << "<td>";
    if (href.size() == 0)
	addImage(file, caption);
    else
	addImage(file, caption, href);
    out_ << "<br/>" << caption << "</td>" << endl;
}

// ---------------------------------------------------------------------------------------
void HtmlReport::addImage(std::string file, std::string caption)
{
    out_ << "<img src=\"" << file << "\" alt=\"" << caption << "\"/>";
}

// ---------------------------------------------------------------------------------------
void HtmlReport::addImage(std::string file, std::string caption, std::string href)
{
    out_ << "<a href=\"" << href << "\"><img src=\"" << file << "\" alt=\"" << caption << "\"/></a>";
}

// ---------------------------------------------------------------------------------------
void HtmlReport::addH(std::string text, char level, bool foldable)
{
    flushDiv();
    if (isTable_) endTable();
    out_ << "<h" << level;
    if (foldable) out_ << " class=\"expanded\"";
    isCollapsable_ = foldable;
    out_ << ">" << text << "</h" << level << ">" << endl;
    return;
}

void HtmlReport::addP(std::string text, bool foldable)
{
    if (isTable_) endTable();
    out_ << "<p";
    if (foldable) out_ << " class=\"expanded\"";
    isCollapsable_ = foldable;
    out_ << ">" << text << "</p>" << endl;
    return;
}

// ---------------------------------------------------------------------------------------
void HtmlReport::beginDiv()
{
    if (isCollapsable_)
	out_ << "<div class=\"show\">" << endl;
    else
	out_ << "<div>" << endl;
    divCounter_++;
}

void HtmlReport::endDiv()
{
    if (divCounter_ > 0)
    {
	out_ << "</div>" << endl;
	divCounter_--;
    }
}

void HtmlReport::flushDiv()
{
    while (divCounter_ > 0)
	endDiv();
}

