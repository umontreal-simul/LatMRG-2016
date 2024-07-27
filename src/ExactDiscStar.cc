#include "latmrg/ExactDiscStar.h"
#include <cassert>
#include <cmath>
#include <string>

using namespace std;


namespace LatMRG
{

//===========================================================================

void ExactDiscStar::init (int d)
{
   assert ((d == 3) || (d == 2));
   m_N = m_n * m_r;
   m_Ksi = new double[m_n + 2];
   m_y = new double[m_n + 2];
   m_x = new double[m_n + 2];
   m_frac = new double[m_n + 2];

   if (3 == d) {
   m_Eta = new double[m_n + 2];
   m_z = new double[m_n + 2];
   m_temp = new double[m_n + 2];
   }
}


//===========================================================================

ExactDiscStar::ExactDiscStar (long n, bool prime, int d) :
      Discrepancy (n, prime, 0, d)
{
   init (d);
}


//===========================================================================

ExactDiscStar::~ExactDiscStar ()
{
   delete[] m_x;
   delete[] m_y;
   delete[] m_Ksi;
   delete[] m_frac;

   if (3 == m_dim) {
   delete[] m_z;
   delete[] m_Eta;
   delete[] m_temp;
   }
}


//===========================================================================

const string ExactDiscStar::toString () const
{
   string S(Discrepancy::toString ());
   string::size_type i = S.find("gamma");
   string::size_type j = S.find("\n", i);
   S.erase(i, j);
   return "ExactDiscStar:\n" + S;
}


//===========================================================================

void ExactDiscStar::insert2Dim (double y, int n)
{
   // Insert y in sorted array Ksi so that Ksi is always sorted
   if (0 == n) {
      m_Ksi[1] = m_Ksi[0];
      m_Ksi[0] = y;
      return;
   }
   int i = n;
   while (y < m_Ksi[i]) {
      m_Ksi[i + 1] = m_Ksi[i];
      i--;
   }
   m_Ksi[i + 1] = y;
}


//===========================================================================

void ExactDiscStar::insert_yz (double y, double z, int n)
{
   // Insert y in sorted array Ksi so that Ksi is always sorted
   int i = n;
   while (y < m_Ksi[i]) {
      m_Ksi[i + 1] = m_Ksi[i];
      m_temp[i + 1] = m_temp[i];
      i--;
   }
   m_Ksi[i + 1] = y;
   m_temp[i + 1] = z;
}


//===========================================================================

void ExactDiscStar::insert_z (double z, int n)
{
   // Insert y in sorted array Ksi so that Ksi is always sorted
   int i = n;
   while (z < m_Eta[i]) {
      m_Eta[i + 1] = m_Eta[i];
      i--;
   }
   m_Eta[i + 1] = z;
}


//===========================================================================

double ExactDiscStar::compute (long z[], int d)
{
   if (d == 2)
      return compute2Dim (z);
   if (d == 3)
      return compute3Dim (z);
   assert ((d == 2) || (3 == d));
   return -1.0e100;
}


//===========================================================================

double ExactDiscStar::compute2Dim (long z[])
{
   long s;
   const long zd = z[2];
   long k;
   long i;

   for (i = 1; i <= m_n; i++) {
      s = ((i - 1) * zd) % m_n;
      m_y[i] = (double) s / m_n;
      m_x[i] = (double) (i - 1) / m_n;
      m_frac[i] = (double) i / m_n;
   }
   m_x[0] = m_y[0] = 0;
   m_x[m_n + 1] = m_y[m_n + 1] = 1;
   m_frac[0] = 0;
   m_frac[m_n] = 1;
   m_Ksi[0] = 1;

   double tem, sem;
   double D = 0;               // discrepancy
   for (i = 0; i <= m_n; i++) {
      insert2Dim (m_y[i], i);
      for (k = 0; k <= i; k++) {
         tem = fabs (m_frac[k] - m_x[i] * m_Ksi[k]);
         sem = fabs (m_frac[k] - m_x[i + 1] * m_Ksi[k + 1]);
         if (sem > tem)
            tem = sem;
         if (tem > D)
            D = tem;
      }
   }

   return D;
}


//===========================================================================

double ExactDiscStar::compute3Dim (long z[])
{
   double D = 0.;              // discrepancy
   double tem, rem;
   long s, t;
   const long zd = z[2];
   const long zd3 = z[3];
   long k, i, l;

   for (i = 1; i <= m_n; i++) {
      s = (i * zd) % m_n;
      t = (i * zd3) % m_n;
      m_y[i] = (double) s / m_n;
      m_z[i] = (double) t / m_n;
      m_x[i] = (double) i / m_n;
   }
   m_x[0] = m_y[0] = m_z[0] = 0.;
   m_x[m_n + 1] = m_y[m_n + 1] = m_z[m_n + 1] = 1.;
   m_Ksi[0] = 0.;
   m_Ksi[1] = 1.;
   m_Eta[0] = 0.;
   m_Eta[1] = 1.;

   for (i = 0; i <= m_n; i++) {
      insert_yz (m_y[i], m_z[i], i);
      for (k = 0; k <= i; k++) {
         insert_z (m_temp[k], k);
         for (l = 0; l <= k; l++) {
            tem = fabs (m_x[l] - m_x[i] * m_Ksi[k] * m_Eta[l]);
            rem = fabs (m_x[l] - m_x[i + 1] * m_Ksi[k + 1] * m_Eta[l + 1]);
            if (rem > tem)
               tem = rem;
            if (tem > D)
               D = tem;
         }
      }
   }

   return D;
}


//===========================================================================

}
