#include "latmrg/PalphaLCG.h"
#include "latcommon/Util.h"
#include "latcommon/Num.h"
#include "latmrg/ParamReader.h"
#include "latmrg/LatConfig.h"

extern "C" {
#include "unif01.h"
#include "ulcg.h"
#include "num.h"
}

#include <cassert>
#include <cmath>


using namespace LatCommon;


namespace
{
const double Pi = 3.14159265358979323844;
const double DeuxPi2 = 2.0 * Pi * Pi;
const double DeuxPi4 = DeuxPi2 * Pi * Pi;
const double DeuxPi6 = DeuxPi4 * Pi * Pi;
const double DeuxPi8 = DeuxPi6 * Pi * Pi;
const double UnSixieme = 1.0 / 6.0;


void trace (char *mess)
{
   std::cout << "*** " << mess << std::endl;
}

}


namespace LatMRG
{

//========================================================================
//Procedure qui lit les caracteristiques du generateur dans le fichier
PalphaLCG::PalphaLCG (const LatConfig & config)
{
   m_genType = config.genType[0];
   conv (m_m, config.comp[0]->getM());
   m_prime = config.primeM;
   conv (m_a, config.comp[0]->a[1]);
//   long m_infa;                   // Lower bound for m_a
//   long m_supa;                   // Upper bound for m_a
   m_calcType = config.calcPalpha;
   m_minDim = config.td[0];
   m_maxDim = config.td[1];
   m_alpha = config.alpha;  
   m_seed = config.seed;
   CreateVect (TabU, m_maxDim);
   CreateVect (Fac, m_maxDim);
   CreateVect (m_Beta, m_maxDim);
   for (int i = 0; i <= m_maxDim; ++i)
      m_Beta[i] = config.Beta[i];          
}


//========================================================================

PalphaLCG::~PalphaLCG ()
{
   DeleteVect (Fac);
   DeleteVect (TabU);
   DeleteVect (m_Beta);
}


//========================================================================

void PalphaLCG::setMult (long a)
{
   m_a = a;
}


//========================================================================

void PalphaLCG::calcPalpha2 (int minDim, int maxDim, double P2[])
// Calcule P_2 en dimensions S1 a S2, en supposant que les poids sont 1.
// Les valeurs sont retournees dans P2[S1..S2].
{
// trace ("calcPalpha2 --3");
   assert (maxDim <= m_maxDim);
   double OldMult, NewMult, Uij;
   int i, s, PosTab = -1;
   long M, A, Seed;
   double Mult[1 + maxDim];
   double Reste[1 + maxDim];
   unif01_Gen *gen;

   // Initialisation
   M = m_m;
   A = m_a;
   Seed = m_seed;
   gen = ulcg_CreateLCG (M, A, 0, Seed);

   for (s = 0; s < maxDim; ++s) {
      Uij = unif01_StripD (gen, 0);
      TabU[s] = 1.0 + DeuxPi2 * (UnSixieme + Uij * (Uij - 1.0));
   }

   const double FactReste = 1.0 + DeuxPi2 * UnSixieme;
   Mult[0] = 1.0;
   Reste[0] = 1.0;
   for (s = 1; s <= maxDim; ++s) {
      Mult[s] = Mult[s - 1] * TabU[maxDim - s];
      P2[s] = Mult[s];
      Reste[s] = Reste[s - 1] * FactReste;
   }

   for (i = 0; i <= m_m - 3; ++i) {
      Uij = unif01_StripD (gen, 0);
      NewMult = 1.0 + DeuxPi2 * (UnSixieme + Uij * (Uij - 1.0));
      for (s = minDim; s <= maxDim; ++s) {
         PosTab = (i - s) % maxDim;
         if (PosTab < 0)
            PosTab += maxDim;
         OldMult = TabU[PosTab];
         Mult[s] *= NewMult / OldMult;
         P2[s] += Mult[s];
      }
      TabU[PosTab] = NewMult;
   }

   for (s = minDim; s <= maxDim; ++s) {
      P2[s] = (P2[s] + Reste[s]) / M - 1.0;
   }

   ulcg_DeleteGen (gen);
}


//========================================================================

double PalphaLCG::calcPalpha2 (int dim)
{
// trace ("calcPalpha2");
   assert (dim <= m_maxDim);
//   assert (1.0 == m_Beta[dim] && 1.0 == m_Beta[1]);
   double NewMult, U;
   int i, j, pos;
   unif01_Gen *gen;

   const long M = m_m;
   const long A = m_a;
   const long Seed = m_seed;
   gen = ulcg_CreateLCG (M, A, 0, Seed);

   double Mult = 1.0;
   double Reste = 1.0;
   for (j = 0; j < dim; ++j) {
      U = unif01_StripD (gen, 0);
      Fac[j] = DeuxPi2;
      TabU[j] = 1.0 + Fac[j] * (UnSixieme + U * (U - 1.0));
      Mult = Mult * TabU[j];
      Reste = Reste * (1.0 + Fac[j] * UnSixieme);
   }
   double Add = Mult;

   for (i = 0; i < m_m - 2; ++i) {
      U = unif01_StripD (gen, 0);
      pos = i % dim;
      NewMult = 1.0 + Fac[pos] * (UnSixieme + U * (U - 1.0));
      Mult = Mult * NewMult / TabU[pos];
      TabU[pos] = NewMult;
      Add = Add + Mult;
   }

   ulcg_DeleteGen (gen);
   Add = (Add + Reste) / M - 1.0;
   return Add;
}


//========================================================================

void PalphaLCG::calcPalpha2PerNonMax (int minDim, int maxDim, double P2[])
// Calcule P_2 en dimensions S1 a S2, en supposant que les poids sont 1.
// Les valeurs sont retournees dans P2[S1..S2].
// Le tableau P2 doit bien sur commencer a 0.
{
// trace ("calcPalpha2PerNonMax --3");
   assert (maxDim <= m_maxDim);
   double Prod, U;
   int j, s;
   long M, A;
   long Ind[1 + maxDim];
   long Puis[1 + maxDim];
   double *T;

   // Initialisation
   M = m_m;
   A = m_a;
   const double InvM = 1.0 / M;
   CreateVect (T, M - 1);
   for (j = 0; j < M; ++j) {
      U = j * InvM;
      T[j] = 1.0 + DeuxPi2 * (UnSixieme + U * (U - 1.0));
   }

   Puis[0] = 1;
   for (s = 1; s <= maxDim; ++s) {
      Puis[s] = num_MultModL (Puis[s - 1], A, 0, M); // = A^s mod M
      Ind[s] = 0;
   }

   for (j = 0; j < M; ++j) {
      Prod = 1.0;
      for (s = 1; s <= maxDim; ++s) {
         Prod *= T[Ind[s]];
         P2[s] += Prod;
         Ind[s] = (Ind[s] + Puis[s]) % M;
//         assert (Ind[s] >= 0);
      }
   }

   for (s = minDim; s <= maxDim; ++s) {
      P2[s] = P2[s] / M - 1.0;
   }
   DeleteVect (T);
}


//========================================================================

double PalphaLCG::calcPalpha2PerNonMax (int dim)
{
// trace ("calcPalpha2PerNonMax");
   assert (dim <= m_maxDim);
   const long M = m_m;
   const long A = m_a;
   int i;
   long Vect[dim + 1];
   long VectTemp[dim + 1];

   // Tableau contenant le vecteur des points
   Vect[0] = 1;
   VectTemp[0] = 1;
   Fac[0] = DeuxPi2 * m_Beta[1] * m_Beta[1];
   double Reste = 1.0 + Fac[0] / 6.0;

   for (i = 1; i < dim; ++i) {
      Vect[i] = num_MultModL (Vect[i - 1], A, 0, M);
      VectTemp[i] = Vect[i];
      Fac[i] = DeuxPi2 * m_Beta[i + 1] * m_Beta[i + 1];
      Reste = Reste * (1.0 + Fac[i] / 6.0);
   }

   // Calcul du P_alpha
   const double MM = M;
   double Mult = 1.0;
   double Add = 0.0;
   double U;
   for (i = 1; i < M; ++i) {
      Mult = 1.0;
      for (int j = 0; j < dim; ++j) {
         U = VectTemp[j] / MM;
         Mult = Mult * (1.0 + Fac[j] * (UnSixieme + U * (U - 1.0)));
      }
      Add += Mult;
      for (int k = 0; k < dim; ++k) {
         VectTemp[k] = (VectTemp[k] + Vect[k]) % M;
      }
   }
   Add = (((Add + Reste) / MM) - 1.0) * m_Beta[0];
   return Add;
}


//========================================================================

double PalphaLCG::calcPalpha4 (int dim)
{
// trace ("calcPalpha4");
   assert (dim <= m_maxDim);
//   assert (1.0 == m_Beta[dim] && 1.0 == m_Beta[1]);
   double U;
   const double UNTRENTE = 1.0 / 30.0;
   const long M = m_m;
   const long A = m_a;
   const long Seed = m_seed;
   unif01_Gen *gen;

   gen = ulcg_CreateLCG (M, A, 0, Seed);

   double Mult = 1.0;
   double Reste = 1.0;
   for (int j = 0; j < dim; ++j) {
      U = unif01_StripD (gen, 0);
      Fac[j] = DeuxPi4 / 3.0;
      TabU[j] = 1.0 - Fac[j] * (((U - 2.0) * U + 1.0)*U*U - UNTRENTE);
      Mult = Mult * TabU[j];
      Reste = Reste * (1.0 + Fac[j] / 30.0); // 1 + (2Pi^4/90) B^4
   }
   double Add = Mult;

   double OldMult;
   int pos;
   for (int i = 0; i < m_m - 2; ++i) {
      U = unif01_StripD (gen, 0);
      pos = i % dim;
      OldMult = TabU[pos];
      TabU[pos] = 1.0 - Fac[pos] * (((U - 2.0) * U + 1.0)*U*U - UNTRENTE);
      Mult = Mult * TabU[pos] / OldMult;
      Add = Add + Mult;
   }

   Add = (Add + Reste) / M - 1.0;
   ulcg_DeleteGen (gen);
   return Add;
}


//========================================================================

double PalphaLCG::calcPalpha4PerNonMax (int dim)
{
// trace ("calcPalpha4PerNonMax");
   assert (dim <= m_maxDim);
   int i;
   const long M = m_m;
   const long A = m_a;
   const double UNTRENTE = 1.0 / 30.0;
   long Vect[dim + 1];
   long VectTemp[dim + 1];

   // Tableau contenant le vecteur des points
   Vect[0] = 1;
   VectTemp[0] = 1;
   Fac[0] = m_Beta[1] * m_Beta[1];
   Fac[0] = DeuxPi4 * Fac[0] * Fac[0] / 3.0;
   double Reste = 1.0 + Fac[0] / 30.0;
   for (i = 1; i < dim; ++i) {
      VectTemp[i] = Vect[i] = num_MultModL (Vect[i - 1], A, 0, M);
      Fac[i] = m_Beta[i + 1] * m_Beta[i + 1];
      Fac[i] = DeuxPi4 * Fac[i] * Fac[i] / 3.0;
      Reste = Reste * (1.0 + Fac[i] / 30.0);
   }

   const double MM = M;
   double Mult = 1.0, Add = 0;
   double U;
   for (i = 1; i < M; ++i) {
      Mult = 1.0;
      for (int j = 0; j < dim; ++j) {
         U = VectTemp[j] / MM;
         Mult *= 1.0 - Fac[j] * (((U - 2.0) * U + 1.0)*U*U - UNTRENTE);
      }
      Add = Add + Mult;
      for (int k = 0; k < dim; ++k) {
         VectTemp[k] = (Vect[k] + VectTemp[k]) % M;
      }
   }

   Add = (((Add + Reste) / M) - 1.0) * m_Beta[0];
   return Add;
}


//========================================================================

double PalphaLCG::calcPalpha6 (int dim)
{
// trace ("calcPalpha6");
   assert (dim <= m_maxDim);
//   assert (1.0 == m_Beta[dim] && 1.0 == m_Beta[1]);
   double U;
   const double QUARAN = 1.0 / 42.0;
   const long M = m_m;
   const long A = m_a;
   const long Seed = m_seed;
   unif01_Gen *gen = ulcg_CreateLCG (M, A, 0, Seed);
   double Mult = 1.0;
   double Reste = 1.0;

   for (int j = 0; j < dim; ++j) {
      Fac[j] = 2.0 * DeuxPi6 / 45.0;
      U = unif01_StripD (gen, 0);
      TabU[j] = 1.0 + Fac[j] *
                ((((U - 3.0) * U + 2.5) * U*U - 0.5) * U*U + QUARAN);
      Mult = Mult * TabU[j];
      Reste = Reste * (1.0 + Fac[j] / 42.0);
   }
   double Add = Mult;

   double OldMult;
   int pos;
   for (int i = 0; i < m_m - 2; ++i) {
      pos = i % dim;
      U = unif01_StripD (gen, 0);
      OldMult = TabU[pos];
      TabU[pos] = 1.0 + Fac[pos] *
            ((((U - 3.0) * U + 2.5) * U*U - 0.5) * U*U + QUARAN);
      Mult = Mult * TabU[pos] / OldMult;
      Add = Add + Mult;
   }
   Add = ((Add + Reste) / M - 1.0);
   ulcg_DeleteGen (gen);
   return Add;
}


//======================================================================

double PalphaLCG::calcPalpha6PerNonMax (int dim)
{
// trace ("calcPalpha6PerNonMax");
   assert (dim <= m_maxDim);
   int i;
   const double QUARAN = 1.0 / 42.0;
   const long M = m_m;
   const long A = m_a;
   long Vect[dim + 1];
   long VectTemp[dim + 1];

   // Tableau contenant le vecteur des points
   Vect[0] = 1;
   VectTemp[0] = 1;
   double z = m_Beta[1] * m_Beta[1];
   z *= m_Beta[1];
   Fac[0] = 2.0 * DeuxPi6 * z * z / 45.0;
   double Reste = 1.0 + Fac[0] * QUARAN;
   for (i = 1; i < dim; ++i) {
      VectTemp[i] = Vect[i] = num_MultModL (Vect[i - 1], A, 0, M);
      z = m_Beta[i + 1] * m_Beta[i + 1];
      Fac[i] = 2.0 * DeuxPi6 * z * z * z / 45.0;
      Reste = Reste * (1.0 + Fac[i] * QUARAN);
   }

   const double MM = M;
   double Mult = 1.0, Add = 0;
   double U;
   for (i = 1; i < M; ++i) {
      Mult = 1.0;
      for (int j = 0; j < dim; ++j) {
         U = VectTemp[j] / MM;
         Mult *= (1.0 + Fac[j] *
                   ((((U - 3.0) * U + 2.5) * U*U - 0.5) * U*U + QUARAN));
      }
      Add = Add + Mult;
      for (int k = 0; k < dim; ++k) {
         VectTemp[k] = (Vect[k] + VectTemp[k]) % M;
      }
   }

   Add = ((Add + Reste) / MM - 1.0) * m_Beta[0];
   return Add;
}


//========================================================================

double PalphaLCG::calcPalpha8 (int dim)
{
// trace ("calcPalpha8");
   assert (dim <= m_maxDim);
//   assert (1.0 == m_Beta[dim] && 1.0 == m_Beta[1]);
   double OldMult;
   double U;
   int pos;
   const double UNTRENTE = 1.0 / 30.0;
   const double DTIERS = 2.0 / 3.0;
   const double STIERS = 7.0 / 3.0;
   const double QTIERS = 14.0 / 3.0;

   const long M = m_m;
   const long A = m_a;
   const long Seed = m_seed;
   unif01_Gen *gen = ulcg_CreateLCG (M, A, 0, Seed);

   double Mult = 1.0;
   double Reste = 1.0;
   for (int j = 0; j < dim; ++j) {
      U = unif01_StripD (gen, 0);
      Fac[j] = DeuxPi8 / 315.0;
      TabU[j] = 1.0 - Fac[j] * (((((U - 4.0) * U + QTIERS) *
                 U*U - STIERS) * U*U + DTIERS) * U*U - UNTRENTE);
      Mult = Mult * TabU[j];
      Reste = Reste * (1.0 + Fac[j] * UNTRENTE);
   }
   double Add = Mult;

   for (int i = 0; i < m_m - 2; ++i) {
      U = unif01_StripD (gen, 0);
      pos = i % dim;
      OldMult = TabU[pos];
      TabU[pos] = 1.0 - Fac[pos] * (((((U - 4.0) * U +
                     QTIERS) * U*U - STIERS) * U*U + DTIERS) * U*U - UNTRENTE);
      Mult = Mult * TabU[pos] / OldMult;
      Add = Add + Mult;
   }
   ulcg_DeleteGen (gen);

   Add = ((Add + Reste) / M - 1.0);
   return Add;
}


//========================================================================

double PalphaLCG::calcPalpha8PerNonMax (int dim)
{
// trace ("calcPalpha8PerNonMax");
   assert (dim <= m_maxDim);
   const long M = m_m;
   const long A = m_a;
   int i;
   const double UNTRENTE = 1.0 / 30.0;
   const double DTIERS = 2.0 / 3.0;
   const double STIERS = 7.0 / 3.0;
   const double QTIERS = 14.0 / 3.0;
   long Vect[dim + 1];
   long VectTemp[dim + 1];

   // Tableau contenant le vecteur des points
   Vect[0] = 1;
   VectTemp[0] = 1;
   double z = m_Beta[1] * m_Beta[1];
   z = z * z;
   Fac[0] = DeuxPi8 * z * z / 315.0;
   double Reste = 1.0 + Fac[0] * UNTRENTE;
   for (i = 1; i < dim; ++i) {
      Vect[i] = num_MultModL (Vect[i - 1], A, 0, M);
      VectTemp[i] = Vect[i];
      z = m_Beta[i + 1] * m_Beta[i + 1];
      z = z * z;
      Fac[i] = DeuxPi8 * z * z / 315.0;
      Reste = Reste * (1.0 + Fac[i] * UNTRENTE);
   }

   // Calcul du P_alpha
   const double MM = M;
   double Mult = 1.0;
   double Add = 0.0;
   double U;
   for (i = 1; i < M; ++i) {
      Mult = 1.0;
      for (int j = 0; j < dim; ++j) {
         U = VectTemp[j] / MM;
         Mult *= (1.0 - Fac[j] * (((((U - 4.0) * U +
                  QTIERS) * U*U - STIERS) * U*U + DTIERS) * U*U - UNTRENTE));
      }
      Add = Add + Mult;
      for (int k = 0; k < dim; ++k) {
         VectTemp[k] = (Vect[k] + VectTemp[k]) % M;
      }
   }

   Add = ((Add + Reste) / MM - 1.0) * m_Beta[0];
   return Add;
}


//====================================================**************

double PalphaLCG::calcPalpha2Verif (int dim)
{
// trace ("calcPalpha2Verif");
   assert (dim <= m_maxDim);
   double Mult;
   double Add;
   int i, j;
   double U;
   int pos;
   long M, A, Seed;
   double Reste;
   unif01_Gen *gen;

   M = m_m;
   A = m_a;
   Seed = m_seed;
   gen = ulcg_CreateLCG (M, A, 0, Seed);

   // Declaration d'un tableau de grandeur dim pour garder
   // les U qui serviront dans le deuxieme for (

   // Initialisation
   Mult = 1.0;
   Reste = 1.0;
   for (j = 0; j < dim; ++j) {
      U = unif01_StripD (gen, 0);
      TabU[j] = U;
      Mult = Mult * (1.0 + DeuxPi2 * (1.0 / 6.0 + U * (U - 1.0)));
      Reste = Reste * (1.0 + DeuxPi2 / 6.0);
   }
   Add = Mult;

   for (i = 0; i < m_m - 2; ++i) {
      pos = i % dim;
      Mult = 1.0;
      for (j = 0; j < dim; ++j) {
         if (j != pos) {
            Mult *= 1.0 + DeuxPi2 * (1.0 / 6.0 + TabU[j] * (TabU[j] - 1.0));
         }
      }
      U = unif01_StripD (gen, 0);

      Mult = Mult * (1.0 + DeuxPi2 * (1.0 / 6.0 + U * (U - 1.0)));
      TabU[pos] = U;
      Add = Add + Mult;
   }

   Add = (Add + Reste) / M - 1.0;
   ulcg_DeleteGen (gen);
   return Add;
}

//========================================================================
}
