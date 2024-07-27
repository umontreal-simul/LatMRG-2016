/*
 DoubleFormatter.cc for ISO C++
 version 1.00
 modified 16/06/06

 Authors: Hicham Wahbi
          Frédérik Rozon
          Richard Simard
*/

#include <sstream>
#include <string>
#include <iomanip>

#include "latmrg/DoubleFormatter.h"

using namespace std;

namespace LatMRG
{
DoubleFormatter::DoubleFormatter(int prec)
{
   m_precision = prec;
}

string DoubleFormatter::format(void* v)
{
   double val = *static_cast<double*>(v);
   ostringstream outs;

   outs << setprecision(m_precision) << val;

   return outs.str();
}
}
