/*
 ReportLat.cc for ISO C++
 version 1.00
 
 authors: Hicham Wahbi
          Frédérik Rozon
          Richard Simard
*/

#include "latmrg/ReportLat.h"
#include "latmrg/TableColumnImpl.h"
#include "latmrg/Writer.h"

using namespace std;
using namespace LatCommon;

namespace LatMRG
{
ReportLat::ReportLat(Writer* writer, LatConfig* config, ReportHeader* header,
                     ReportFooter* footer): m_dFormat(6)

{
   m_config = config;
   m_writer = writer;
   m_header = header;
   m_footer = footer;

   TableColumnImpl<int>* dims = new TableColumnImpl<int>(&m_iFormat, "t");

   m_results.add(dims);

   for (int i = m_config->td[0]; i <= m_config->td[1]; i++) {
      dims->addValue(static_cast<void*>(&i));
   }

   m_base_col = 1;
}

ReportLat::~ReportLat()
{
}

void ReportLat::printHeader()
{
   m_header->printHeader();
}

void ReportLat::printFooter()
{
   m_footer->printFooter();
}

void ReportLat::printTable()
{
   m_writer->writeTable(m_results, "llllllllll");
}

void ReportLat::baseUpdate(Base & base)
{
   m_writer->writeString (base.toString());
}

void ReportLat::baseUpdate(Base & base, int i)
{
   m_writer->writeString (base.toString(i));
}

void ReportLat::resultUpdate(double results[], int n)
{
   for (int i = 0; i < n; i++) {
      m_results[m_base_col + i]->addValue(static_cast<void*>(&results[i]));
   }
}

void ReportLat::testInit(const string & test, string headers[], int n)
{
   TableColumnImpl<double>* col;

   for (int i = 0; i < n; i++) {
      col = new TableColumnImpl<double>(&m_dFormat, headers[i]);
      m_results.add(col);
   }

   m_next_base = m_base_col + n;

//   m_writer->writeString("Test: ");
//   m_writer->writeString(test);
//   m_writer->newLine();
}

void ReportLat::testCompleted()
{
   m_writer->newLine();
   m_writer->writeString("Test completed successfully");
   m_writer->newLine();

   m_base_col = m_next_base;
}

void ReportLat::testFailed(int dim)
{
   m_writer->newLine();
   m_writer->writeString("Test failed at dim: ");
   m_writer->writeInt(dim);
   m_writer->newLine();

   m_base_col = m_next_base;
}
}
