#include "latmrg/IntFactorization.h"
#include "NTL/ZZ.h"
#include "latcommon/Util.h"

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <cassert>

using namespace LatCommon;

#define USE_MIRACL

namespace LatMRG
{

//===========================================================================

IntFactorization::IntFactorization (const char *name):
      m_status(UNKNOWN)
{
   if (name != 0)
      read(name);
   else
      m_number = 0;
}


//===========================================================================

IntFactorization::IntFactorization (const MScal & n):
      m_number (n), m_status(UNKNOWN)
{}


//===========================================================================

IntFactorization::IntFactorization (const IntFactorization & f)
      : m_number(f.m_number), m_status(f.m_status), m_factorList(f.m_factorList),
      m_invFactorList(f.m_invFactorList)
{}


//===========================================================================

IntFactorization & IntFactorization::operator= (
   const IntFactorization & f)
{
   if (this != &f) {
      m_factorList = f.m_factorList;
      m_invFactorList = f.m_invFactorList;
      m_number = f.m_number;
      m_status = f.m_status;
   }
   return *this;
}


//===========================================================================

void IntFactorization::clear ()
{
   m_status = UNKNOWN;
   m_number = 0;
   m_factorList.clear ();
   m_invFactorList.clear ();
}

//===========================================================================

IntFactorization::~IntFactorization ()
{
   clear ();
}


//===========================================================================

void IntFactorization::addFactor (const MScal & x, int mult,
                                  PrimeType st)
{
   IntFactor f (x, mult, st);
   m_factorList.push_back (f);
}


//===========================================================================

void IntFactorization::read (const char *name) throw (std::invalid_argument)
{
   std::ifstream in (name);
   if (!in) {
      std::string str("IntFactorization::read:   Unable to open input file  ");
      str += name;
      throw std::invalid_argument(str);
   }

   in >> m_number;   // read the m_number
   in.ignore (100, '\n'); // drop rest of line
   
   int vsize = 0;
   std::string tampon;

   while (in >> tampon) {
      MScal x;
      int k;
      char c;
      PrimeType stat;
      NTL::ZZ f;

      f = NTL::to_ZZ (tampon.c_str ());
      conv (x, f);
      in >> k;
      in >> c;

      switch (c) {
      case 'P':
         stat = PRIME;
         break;
      case 'Q':
         stat = PROB_PRIME;
         break;
      case 'C':
         stat = COMPOSITE;
         break;
      default:
         stat = UNKNOWN;
      }

      addFactor (x, k, stat);
      ++vsize;
   }

   //unique ();
   assert (true == checkProduct ());
     
   m_invFactorList.reserve (vsize);
}

//===========================================================================

void IntFactorization::unique ()
{
   sort ();
   int j = 1;

   std::list<IntFactor>::iterator it = m_factorList.begin ();
   std::list<IntFactor>::iterator it2 = m_factorList.begin ();
   if (it2 != m_factorList.end ())
      ++it2;
   while (it2 != m_factorList.end ()) {
      if (it->getFactor () == it2->getFactor ()) {
         it->setMultiplicity (++j);
         it2 = m_factorList.erase (it2);
      } else {
         j = 1;
         ++it;
         ++it2;
      }
   }
}


//===========================================================================

std::string IntFactorization::toString () const
{
   std::list<IntFactor>::const_iterator it =
      m_factorList.begin ();
   std::ostringstream out;
   out << m_number << std::endl;
//   out << "its factors:\n";
   while (it != m_factorList.end ()) {
      out << "\t" << (*it).toString () << std::endl;
      ++it;
   }
   out << std::endl;
/*
   if (!m_invFactorList.empty()) {
      out << "the inverse factors:\n";
      for (unsigned int i = 0; i < m_invFactorList.size(); ++i) {
         out << m_invFactorList[i] << std::endl;
      }
   }
*/
   return out.str ();
}


//===========================================================================

void IntFactorization::factorize ()
{
#ifndef USE_MIRACL
   NTL::Error (
      "IntFactorization::factorize -->   gdef::USE_MIRACL undefined ");
#else

   char *prog;          // Home of factorizing program
   prog = getenv ("FACTORHOME");
   if (!prog) {
      NTL::Error (
         "IntFactorization::factorize -->   FACTORHOME is not defined");
   }
   std::string S(prog);
   S += "/factorm ";
   std::ostringstream num;
   num << m_number;
   S += num.str();

   // Choose a temporary name for the file
   const char *filename = "temp938573291";

#if 0

   std::ofstream file (filename);
   std::streambuf *buffer = std::cout.rdbuf ();
   // Redirect standard output to file
   std::cout.rdbuf (file.rdbuf ());
   // factorize; this may take a very long time
   std::system (S.c_str ());
   // Redirect standard output to cout
   std::cout.rdbuf (buffer);
#endif

   S += " > ";
   S += filename;
   // factorize and set output to filename
   std::system (S.c_str ());

   // Now read the result file and extract the prime factors from the
   // lines PRIME FACTOR xxx
   std::ifstream in (filename);
   if (!(in.is_open())) {
      std::cerr << "Error:   cannot open file   filename\n";
      exit(8);
   }
   std::string line;
   std::string::size_type pos;
   MScal z;
   getline (in, line, '\n');
   pos = line.find ("this number is prime");
   if (pos != std::string::npos) {
      addFactor (m_number, 1, PRIME);
      remove
         (filename);
      return;
   }
   while (getline (in, line, '\n')) {
      pos = line.find ("PRIME FACTOR");
      if (pos != std::string::npos) {
         // Found a prime factor
         S = line.substr (pos + 12);
         conv(z, S.c_str ());
         addFactor (z, 1, PRIME);
      }
      pos = line.find ("COMPOSITE FACTOR");
      if (pos != std::string::npos) {
         // The program gave up: last factor is composite
         S = line.substr (pos + 16);
         conv(z, S.c_str ());
         addFactor(z, 1, COMPOSITE);
      }
      pos = line.find ("MIRACL error");
      if (pos != std::string::npos) {
         // error in factorizing
         Error ("MIRACL error in IntFactorization::factorize:   Number too big");
      }
   }
   unique ();
   remove (filename);
#endif

}


//===========================================================================

bool IntFactorization::checkProduct () const
{
   MScal temp;
   temp = 1;

   std::list<IntFactor >::const_iterator it =
      m_factorList.begin ();

   while (it != m_factorList.end ()) {
      for (int j = it->getMultiplicity (); j > 0; --j)
         temp *= it->getFactor ();
      ++it;
   }
   if (temp != m_number) {
      std::cout << "checkProduct:   ERROR -->  " << m_number
      << " != " << temp << std::endl;
   }
   return (temp == m_number);
}


//===========================================================================

void IntFactorization::calcInvFactors ()
{
   //   sort ();
   unsigned int j = 0;
   std::list<IntFactor>::reverse_iterator it =
      m_factorList.rbegin ();

   while (it != m_factorList.rend ()) {
      ++j;
      if (m_invFactorList.capacity () < j) {
         m_invFactorList.clear ();
         m_invFactorList.reserve (j + 10);
      }
      m_invFactorList.push_back (m_number / it->getFactor ());
      ++it;
   }
}

//===========================================================================

}                                 // namespace LatMRG
