/*
 ReportFooterLat.cc for ISO C++
 version 1.00
 authors: Hicham Wahbi
          Frédérik Rozon
          Richard Simard
*/

#include "latmrg/Merit.h"
#include "latmrg/ReportFooterLat.h"
#include "latcommon/Const.h"


namespace LatMRG
{
ReportFooterLat::ReportFooterLat (Writer * writer, LatticeTest * test):
      ReportFooter(writer)
{
   m_test = test;
}


void ReportFooterLat::printFooter()
{
   int dimMin = m_test->getMinDim ();
   int dimMax = m_test->getMaxDim ();
   int dimWorst;
   // Get minimal merit value and dimension where it occurs
   double S = m_test->getMerit().getST (dimMin, dimMax, dimWorst);
   if (m_test->getLattice()->getNorm() == LatCommon::L2NORM)
      S = sqrt(S);
   m_writer->newLine ();
   m_writer->writeString (" Min merit:   S_");
   m_writer->writeInt (dimWorst);
   m_writer->writeString (" = ");
   m_writer->writeDouble (S);
   m_writer->newLine ();
   m_writer->newLine ();
}

}
