#include <fstream>
#include <iostream>
#include <vector>

#include "latcommon/Types.h"
#include "latcommon/Util.h"
#include "latcommon/Const.h"
#include "latmrg/Chrono.h"

#include "latmrg/SeekConfig.h"
#include "latmrg/ParamReaderSeek.h"
#include "latmrg/IntPrimitivity.h"

#include "latmrg/MRGLattice.h"
#include "latmrg/KorobovLattice.h"
#include "latmrg/Rank1Lattice.h"
#include "latmrg/MRGLatticeFactory.h"
#include "latmrg/MRGComponent.h"

#include "latcommon/NormaBestLat.h"
#include "latcommon/NormaLaminated.h"
#include "latcommon/NormaRogers.h"
#include "latcommon/NormaMinkL1.h"
#include "latcommon/NormaMinkowski.h"

#include "latmrg/LatTestBeyer.h"
#include "latmrg/LatTestSpectral.h"
#include "latmrg/LatTestPalpha.h"
#include "latmrg/ReportHeader.h"
#include "latmrg/Merit.h"
#include "latmrg/WriterRes.h"
#include "latmrg/TestProjections.h"

#include "latmrg/Zone.h"

using namespace std;
// using namespace NTL;
using namespace LatMRG;
using namespace LatCommon;


namespace
{

//===========================================================================

Normalizer *normal;
// Flag to ensure that Normalizer are created only once;
// must know m and k for that.
bool isFirstTest = true;
NormType Norm = L2NORM;

SeekConfig config;
ofstream fout;
// streambuf *psbuf;
Chrono timer;

MRGComponent **compJ;
MMat coef;
LatMRG::IntLattice *lattice;
LatMRG::IntLattice *master;
LatticeTest *latTest;
LatticeTest **pool;            // vecteur de LatticeTest*
int poolLen = 0;               // longueur du vecteur pool <= numGen

Zone **zone;
Zone **reg;
IntPrimitivity **primJ;

long numATried = 0;
// long numAPrimitive = 0;
long *CoNoElemPrim;
long *CoElemPrim;
long *TotEP;
long *CoPerMax;
long *CoRegions;
long *Hk;
long *H;
double *minVal;
// long ExceedBB = 0;    // Number of cases for which m_countNodes > maxNodesBB

void SeekGen (int);



//===========================================================================

void SortBestGen ()
{
   for (int r = 0; r < poolLen - 1; r++) {
      LatticeTest *latTest1 = pool[r];
      double merit1 = latTest1->getMerit ().getWorstMerit ();
      for (int s = r + 1; s < poolLen; s++) {
         LatticeTest *latTest2 = pool[s];
         double merit2 = latTest2->getMerit ().getWorstMerit ();
         if (merit2 > merit1) {
            latTest = pool[s];
            pool[s] = pool[r];
            pool[r] = latTest;
            merit1 = merit2;
         }
      }
   }
}


//===========================================================================

int readConfigFile (int argc, char **argv)
{
   if (argc != 2) {
      cerr << "\n***  Usage: " << argv[0] << " <data file>" << endl;
      return -1;
   }
   // Lecture des paramètres
   string fname (argv[1]);
   fname += ".dat";
   ParamReaderSeek paramRdr (fname.c_str ());
   paramRdr.read (config);

   // Writer *rw;
   switch (config.outputType) {
   case RES:
      fname = argv[1];
      fname += ".res";
      // rw = new WriterRes (fname.c_str());
      fout.open (fname.c_str ());
      break;
   case TEX:
      fname = argv[1];
      fname += ".tex";
      return -2;
      // rw = new WriterTex(fname.c_str());
      break;
   case TERMINAL:
      cerr << "\n*** SeekMain::readConfigFile::outputType: " <<
      "TERMINAL is not implemented." << endl;
      throw std::invalid_argument ("SeekMain::readConfigFile");
      // rw = new WriterRes (&cout);
      // psbuf = cout.rdbuf(); // get cout streambuf
      // fout.rdbuf(psbuf); // assign streambuf to fout
      break;
   default:
      cerr << "\n*** outputType:   no such case" << endl;
      return -2;
   }
   fname.clear ();
   return 0;
}


//===========================================================================

void debug (int dim1, int dim2)
{
   bool racF = false;
   bool invertF = true;
   cout <<
   "=================================================================\n\n";
   for (int j = 0; j < poolLen; j++) {
      LatMRG::IntLattice *lat = pool[j]->getLattice ();
      cout << "A = " << lat->toStringCoef ();
      int dimWorst;
      double S_T = pool[j]->getMerit ().getST (dim1, dim2, dimWorst);
      if (lat->getNorm () == L2NORM)
         racF = true;
      cout << "\t\tS_" << dimWorst << " = " << S_T << endl << endl;
      cout << "   dim\t\t\tl_dim\t\t     S_dim" << endl;
      cout << pool[j]->getMerit ().toString (dim1, dim2, racF,
                                             invertF) << endl;
   }
}


//===========================================================================

void printPool (char *mess)
{
   fout << "--------------------------------------------------" <<
   mess << endl;
   for (int i = 0; i < poolLen; i++) {
      LatticeTest *latTest = pool[i];
      LatMRG::IntLattice *lat = latTest->getLattice ();
      int dimWorst;
      double S_T = latTest->getMerit ().getST (config.td[0],
                   config.td[1], dimWorst);
      if (lat->getNorm () == L2NORM) {
         S_T = sqrt (S_T);
      }
      if (config.J > 1) {
         for (int j = 0; j < config.J; j++) {
            fout << " " << toString (lat->comp[j]->a, 1, lat->comp[j]->k);
            if (j < config.J - 1)
               fout << endl;
         }
      } else {
         fout << " " << lat->toStringCoef ();
      }
      fout << "\t\tS_" << dimWorst << " = " << S_T << endl;
   }
}


//===========================================================================

void PrintComponentData (const Component & comp, int j)
{
   // Prints parameter of component j of combined generator
   fout << "Component " << j + 1 << endl;
   fout << "-----------\n";
   fout << "   GenType                : " << toStringGen (comp.
         genType) << endl;
   fout << "   Modulus m              : " << comp.modulus.m << endl;
   fout << "   Order k                : " << comp.k << endl;
   // fout << " n_j : " << comp.n << endl;
   // fout << " rho_j = (m_j)^k_j-1 : " << comp.rho << endl;
   if (comp.PerMax) {
      fout << "   Factors of m-1         : " <<
      compJ[j]->ifm1.toString ();
      fout << "   Factors of r           : " <<
      compJ[j]->ifr.toString ();
   }
   fout << "   Search method                      : " <<
   toStringSearchMethod (comp.searchMethod) << endl;
   string margin = "   Bounds : ";
   for (int i = 1; i <= comp.k; i++) {
      if (i != 1)
         margin = "            ";
      fout << margin << "a" << i << " from : " << comp.b[i] << endl;
      fout << "                 to : " << comp.c[i] << endl;
   }
   fout << "   Implementation condition           : " <<
   toStringImplemCond (comp.implemCond) << endl;
   fout << "   Maximum period required            : ";
   if (comp.PerMax)
      fout << "yes";
   else
      fout << "no";
   fout << endl << "\n-----------------------------------------------" <<
   endl;
}


//===========================================================================

void PrintPolyStats (const Component & comp, int j)
{
   fout << "Component " << j + 1 << endl;
   fout << "-----------\n";
   if (comp.PerMax) {
      fout << "    Values of a_" << comp.k << " tried                   : "
           << CoNoElemPrim[j] + CoElemPrim[j] << endl;
      fout << "    Values of a_" << comp.k <<
      " primitive element       : " << CoElemPrim[j] << endl;
      fout << "    Polynomials with a_" << comp.k <<
      " primitive element: " << TotEP[j] << endl;
      fout << "    Primitive polynomials                 : "
      << CoPerMax[j] << endl;
   } else
      fout << "    Num. of polynomials examined          : " <<
        TotEP[j] << endl;
   fout << "\n-----------------------------------------------" << endl;
}


//===========================================================================

void PrintResults ()
{
   Component & comp0 = config.compon[0];

   fout << "SEARCH for good ";
   if (comp0.genType == KOROBOV)
      fout << "KOROBOVs" << endl;
   else if (comp0.genType == RANK1)
      fout << "RANK1 Lattices of order " << comp0.k << endl;
   else
      fout << "MRGs of order " << comp0.k << endl;

   // fout << "-----------------------------------------------\n";
   fout << "\nDATA \n";
   fout << "-----------------------------------------------\n";
   for (int j = 0; j < config.J; j++)
      PrintComponentData (config.compon[j], j);

   fout << "   Test                               : "
   << toStringCriterion (config.criter) << endl;
   fout << "   Normalizer                         : "
   << toStringNorma (config.normaType) << endl;
   fout << "   num categories C                   : " << config.C << endl;
   fout << "   min merit                          : " <<
   toString < double *>(config.minMerit, 0, config.C - 1) << endl;
   fout << "   max merit                          : " <<
   toString < double *>(config.maxMerit, 0, config.C - 1) << endl;
   fout << "   d                                  : " << config.d << endl;
   fout << "   td                                 : " <<
   toString (config.td, 0, config.d) << endl;

   if (RANDOM == comp0.searchMethod) {
      fout << "   Seed for RNG                       : " <<
      config.seed << endl;
   }
   fout << "   Max nodes in branch-and-bound      : " <<
   config.maxNodesBB << endl;
   fout << "   Lattice Type                       : " <<
   toStringLattice (config.latType) << endl;
   fout << "   Lattice                            : ";
   if (config.dualF)
      fout << "DUAL" << endl;
   else
      fout << "PRIMAL" << endl;
   fout << "\n\nRESULTS" << endl;
   fout << "-----------------------------------------------\n";
   for (int j = 0; j < -config.J; j++)
      PrintPolyStats (config.compon[j], j);

   fout << "Num. of generators tested             : " <<
          numATried << endl;
   fout << "Num. of generators conserved          : " <<
           poolLen << endl;
   fout << "Total CPU time (after setup)          : " <<
           timer.toString () << endl;
   fout << setprecision (5) << endl << endl;

   bool rac = false;
   if (Norm == L2NORM)
      rac = true;

   for (int i = 0; i < config.C; i++) {
      fout <<
      "+--------------------------------------------------------------------"
      << endl;
      fout << "   " << poolLen << " generators retained for criterion M_";
      for (int j = 1; j < config.d; j++)
         fout << config.td[j] << ",";
      fout << config.td[config.d] << endl;
      fout <<
      "+--------------------------------------------------------------------"
      << endl;
      fout << "   A\t\t\tM_*" << endl;
      fout <<
      "+--------------------------------------------------------------------"
      << endl;
      for (int s = 0; s < poolLen; s++) {
         LatticeTest *latTest = pool[s];
         LatMRG::IntLattice *lat = latTest->getLattice ();
         int dimWorst = latTest->getMerit ().getDimWorst ();
         double S_T = latTest->getMerit ().getWorstMerit ();
         if (rac)
            S_T = sqrt (S_T);
         if (config.J > 1) {
            for (int j = 0; j < config.J; j++) {
               fout << " " << toString (lat->comp[j]->a, 1,
                                        lat->comp[j]->k);
               if (j < config.J - 1)
                  fout << endl;
            }
         } else {
            fout << " " << lat->toStringCoef ();
         }
         fout << "\n S_";
         if (1 == config.d)
            fout << dimWorst;
         else
            fout << "*";
         fout << " = " << S_T << "\n" << endl;
         if (1 == config.d) {
            fout << "   dim\t\t\tl_dim\t\t     S_dim" << endl;
            fout << latTest->getMerit ().toString (config.td[0],
                           config.td[1], rac, config.invertF);
         }
         /* if (config.dualF) fout << lat->getDualBasis ().toString(); else
            fout << "m*" << lat->getPrimalBasis ().toString(); */
         fout <<
         "\n+-----------------------------------------------------------------"
         << endl;
      }
   }
}


//===========================================================================

void Test ()
{
   bool stationary = true;
   Component & comp0 = config.compon[0];

   if (config.J > 1) {         // On doit faire la combinaison des MRG
      for (int s = 0; s < config.J; s++)
         compJ[s]->setA (coef[s]);
      lattice = MRGLatticeFactory::fromCombMRG (compJ, config.J,
                config.getMaxDim (), 0, config.latType, Norm);

   } else {
      if (comp0.genType == MRG || comp0.genType == LCG)
         lattice =
            new MRGLattice (comp0.modulus.m, coef[0],
                            config.getMaxDim (), comp0.k, config.latType, Norm);
      else if (comp0.genType == KOROBOV)
         lattice =
            new KorobovLattice (comp0.modulus.m, coef[0][1],
                                config.getMaxDim (), config.latType, Norm);
      else if (comp0.genType == RANK1) {
         stationary = false;
         lattice = new LatMRG::Rank1Lattice (comp0.modulus.m, coef[0],
                                     config.getMaxDim (), Norm);
      }
   }

   if (isFirstTest) {
      normal = lattice->getNormalizer (config.normaType, config.alpha);
      isFirstTest = false;
   }
   lattice->buildBasis (config.td[0]);

   switch (config.criter) {
   case SPECTRAL:
      latTest = new LatTestSpectral (normal, lattice);
      latTest->setDualFlag (config.dualF);
      latTest->setInvertFlag (config.invertF);
      latTest->setMaxAllDimFlag (true);
      if (1 == config.d) {
         latTest->test (config.td[0], config.td[1], minVal);
         latTest->getMerit ().getST (config.td[0], config.td[1]);
      } else {
         if (comp0.genType == MRG || comp0.genType == LCG)
            master = new MRGLattice (*(MRGLattice *) lattice);
         else if (comp0.genType == KOROBOV)
            master = new KorobovLattice (*(KorobovLattice *) lattice);
         else if (comp0.genType == RANK1)
            master = new LatMRG::Rank1Lattice (*(LatMRG::Rank1Lattice *) lattice);

         master->buildBasis (config.td[1]);
         TestProjections proj (master, lattice, latTest, config.td, config.d);
         proj.setDualFlag (config.dualF);
         proj.setPrintF (false);
         double merit = proj.run (stationary, false, minVal);
         latTest->getMerit ().setWorstMerit (merit);
         delete master;
      }
      break;

   case BEYER:
      latTest = new LatTestBeyer (lattice);
      latTest->setDualFlag (config.dualF);
      latTest->setInvertFlag (config.invertF);
      latTest->setMaxAllDimFlag (true);
      latTest->test (config.td[0], config.td[1], minVal);
      break;

   case PALPHA:
      latTest = new LatTestPalpha (normal, lattice);
      latTest->setDualFlag (config.dualF);
      latTest->setMaxAllDimFlag (true);
      latTest->test (config.td[0], config.td[1], minVal);
      break;

   default:
      fout << "Invalid CriterionType " << config.criter << endl;
      exit (1);
      break;
   }

}


//===========================================================================
#if 0

void VerifyCategories (BVect & Me2)
{
   // Called by InsideExam.  Me2 is either Q2 or S2.
   int MaxI;
   double ToBeat;
   int Ind = Order + 1;
   int IndMin = Ind;
   for (int C = 0; C < NbCat; ++C) {
      ++NumTried[C];
      if ((C + 1 < NbCat) || ((C + 1 == NbCat) && TestCompleted)) {
         // On trouve ou le min est atteint pour cette categorie.
         if ((Criterion == Spectral) && (TabDim[C] <= MaxDimS2))
            MaxI = TabDim[C];
         else
            MaxI = MaxDimS2;
         while (Ind <= MaxI) {
            if (Me2[Ind] <= 0.0)
               return ;
            if (Me2[Ind] < Me2[IndMin])
               IndMin = Ind;
            ++Ind;
         }
         if (Me2[IndMin] > MaxMerit[C] * MaxMerit[C])
            return ;
         WITH PireGen[C] ^ DO
         if ((DimMerit == 0) || ((MaximCat[C]) && (Me2[IndMin] >= Merit))
               || ((!MaximCat[C]) && (Me2[IndMin] <= Merit)))
            ConserverGen (C, IndMin);

      }
   }
}
#endif

//===========================================================================

double FindMinMerit ()
/* Find the minimum merit amongst the gen in category 0, i.e. in pool.
   Returns this merit.
   Side effect: swap that gen with the gen in element 0, so that on return,
   the worst gen is in pool[0] */
{
   double minMerit = 1.0e300;
   int minj = -1;
   for (int j = 0; j < poolLen; j++) {
      double curMerit = pool[j]->getMerit ().getWorstMerit ();
      if (curMerit < minMerit) {
         minj = j;
         minMerit = curMerit;
      }
   }

   if (minj > 0) {
      LatticeTest *t = pool[minj];
      pool[minj] = pool[0];
      pool[0] = t;
   }

   return minMerit;
}


//===========================================================================

void CompareMerit ()
{
   int s;
   if (poolLen < config.numGen[0]) {
      pool[poolLen] = latTest;
      poolLen++;
      FindMinMerit ();

   } else {
      double minMerit = pool[0]->getMerit ().getWorstMerit ();
      double curMerit = latTest->getMerit ().getWorstMerit ();

      if (minMerit < curMerit) {
         LatMRG::IntLattice *lat = pool[0]->getLattice ();
         delete lat;
         delete pool[0];
         pool[0] = latTest;
         lat = latTest->getLattice ();
         minMerit = FindMinMerit ();
         int dim1 = config.td[0];
         int dim2 = config.td[1];
         for (s = dim1; s <= dim2; s++)
            minVal[s] = minMerit;

      } else {
         LatMRG::IntLattice *lat = latTest->getLattice ();
         delete lat;
         delete latTest;
      }
      // printPool (" ESPION_4");
   }
}


//===========================================================================

void ExamThisaj (int j, int i, bool Pow2, ProcII Exam)
{
   // Called by InsideExam.  Examines the current a_j
   Component & comp = config.compon[j];
   if ((i == comp.k) && (coef[j][i] == 0))
      return ;

   bool Ok = true;
   if (!Pow2 && comp.PerMax && (i == comp.k)) {
      // On verifie cond.(1) pour periode max
      if (comp.modulus.primeF) { // m_j is a prime
         if (primJ[j]->isPrimitiveElement (coef[j], comp.k))
            ++CoElemPrim[j];
         else {
            Ok = false;
            ++CoNoElemPrim[j];
         }
      } else {
         Ok = comp.modulus.perMaxPowPrime (coef[j][i]);
         if (Ok)
            ++CoElemPrim[j];
         else
            ++CoNoElemPrim[j];
      }
   }

   if (!Ok)
      return ;

   if (comp.implemCond == EQUAL_COEF) {
      int r = 1;
      while ((r < comp.ncoef) && (i > comp.Icoef[r]))
         ++r;
      r--;
      int s = i - 1;
      while (s > comp.Icoef[r]) {
         coef[j][s] = coef[j][i];
         s--;
      }
      i = s + 1;

   } else if (comp.implemCond == ZERO_COEF) {
      int r = comp.ncoef;
      while ((r > 0) && (i <= comp.Icoef[r]))
         --r;
      i = comp.Icoef[r] + 1;
      if (i < 1)
         i = 1;
   }

   if (i == 1) {
      ++TotEP[j];
      if (Pow2 || (!comp.PerMax) || (comp.k == 1) ||
            compJ[j]->maxPeriod23 (coef[j])) {
         // Note: MaxPeriod verifie cond. (2 et 3) pour per. max.
         // If power of prime, then j = i = Orderj [j] = 1.
         ++CoPerMax[j];
         SeekGen (j + 1);
      }

   } else {
      Exam (j, i - 1);
   }
}


//===========================================================================

void ExamBits (MScal q, int j, int i, int b0, int b1, int NbBits,
               bool Pow2, ProcII Exam)
{
   /* Called by InsideExam (and calls itself recursively) to examine all
      bit patterns with less than NbMaxBits [j] in total, that can be
      obtained by switching from 0 to 1 some of the bits of q from bit b0
      to bit b1. */
   Component & comp = config.compon[j];

   for (int b = b0; b <= b1; b++) {
      coef[j][i] = q + TWO_EXP[b];
      ExamThisaj (j, i, Pow2, Exam);
      if ((b > b0) && (NbBits < comp.NumBits))
         ExamBits (coef[j][i], j, i, b0, b - 1, NbBits + 1, Pow2, Exam);
      coef[j][i] = coef[j][i] - TWO_EXP[b];
      coef[j][i] = coef[j][i] - TWO_EXP[b];
      ExamThisaj (j, i, Pow2, Exam);
      if ((b > b0) && (NbBits < comp.NumBits))
         ExamBits (coef[j][i], j, i, b0, b - 1, NbBits + 1, Pow2, Exam);
      coef[j][i] = coef[j][i] + TWO_EXP[b];
   }
}


//===========================================================================

void InsideExam (Zone * Z, int j, int i, ProcII exam)
{
   /*
      Examines all the possible values for coefficient $a_i$ of component
      \texttt{comp} in this region according to the chosen criteria, and
      calls different methods depending on the criteria. \texttt{exam} is
      the (recursive) method that called \texttt{InsideExam}. */
   Component & comp = config.compon[j];
   MScal q, Temp;
   bool Pow2 = false;
   MScal Eight;
   Eight = 8;
   q = Z->getInf ();

   if (comp.PerMax && (!comp.modulus.primeF) && (comp.modulus.b == 2)
         && comp.modulus.c == 0) {
      Pow2 = true;             // m is a power of 2: want a mod 8 = 5 only.
      Modulo (q, Eight, Temp);
      while (Temp != 5) {
         ++q;
         Modulo (q, Eight, Temp);
      }
   }

   if (comp.implemCond == POWER_TWO) {
      Temp = 0;
      ExamBits (Temp, j, i, 0, comp.HighestBit, 1, Pow2, exam);
      return ;
   }

   Zone::ZoneType No = Z->getNo ();
   MScal sup;
   sup = Z->getSup ();
   // if (q == 1)
   // q = 2;
   while (q <= sup) {
      if (timer.timeOver (config.duration))
         return ;
      if ((comp.implemCond == APP_FACT) && Z->DivQ[No])
         Quotient (comp.modulus.m, q, coef[j][i]);
      else
         coef[j][i] = q;
      ExamThisaj (j, i, Pow2, exam);
      if (Pow2)
         q += Eight;
      else
         ++q;
   }
}


//===========================================================================

void ChoisirBornes (int j)
{
   // Choisir une region au hasard et l'examiner au complet
   Zone *Z = zone[j];
   Zone *R = reg[j];
   ++CoRegions[j];
   R->chooseBoundaries (config.compon[j], Z);
   // cout << (R + 1)->toString();
}


//===========================================================================

void ExamRegion (int j, int i)
{
   // Used in the "Random search" case
   if (timer.timeOver (config.duration))
      return ;
   Zone *R;
   R = i + reg[j];             // On va examiner toute cette region
   InsideExam (R, j, i, ExamRegion);
}


//===========================================================================

void ExamAllZones (int j, int i)
{
   // This method is used in the {\em exhaustive search} case.
   if (timer.timeOver (config.duration))
      return ;
   Zone *Z;
   Z = i + zone[j];
   while (Z != 0) {
      // On va examiner toute cette zone.
      InsideExam (Z, j, i, ExamAllZones);
      Z = Z->nextZone;
   }
}


//===========================================================================

void Init ()
{
   if (MINKL1 == config.normaType)
      Norm = L1NORM;

   pool = new LatticeTest * [config.numGen[0]];
   // memset (pool, 0, config.numGen[0] * sizeof (LatticeTest *));

   coef.SetDims (config.J, 1 + config.getMaxK ());

   CoNoElemPrim = new long[config.J];
   CoElemPrim = new long[config.J];
   TotEP = new long[config.J];
   CoPerMax = new long[config.J];
   CoRegions = new long[config.J];
   Hk = new long[config.J];
   H = new long[config.J];
   minVal = new double[1 + config.getMaxDim ()];
   for (int i = 0; i <= config.getMaxDim (); i++)
      minVal[i] = 0.0;

   compJ = new MRGComponent * [config.J];
   primJ = new IntPrimitivity * [config.J];

   int s;
   for (s = 0; s < config.J; s++) {
      Component & comps = config.compon[s];
      CoNoElemPrim[s] = 0;
      CoElemPrim[s] = 0;
      TotEP[s] = 0;
      CoPerMax[s] = 0;
      CoRegions[s] = 0;

      if (comps.PerMax) {
         compJ[s] =
            new MRGComponent (comps.modulus, comps.k,
                              comps.F1, comps.file1.c_str (),
                              comps.F2, comps.file2.c_str ());
         primJ[s] =
            new IntPrimitivity (compJ[s]->ifm1, comps.modulus.m);

      } else
         compJ[s] = new MRGComponent (comps.modulus.m, coef[s],  comps.k);
   }

   zone = new Zone * [config.J];
   for (s = 0; s < config.J; s++)
      zone[s] = new Zone[1 + config.compon[s].k];

   reg = new Zone * [config.J];
   for (s = 0; s < config.J; s++)
      reg[s] = new Zone[1 + config.getMaxDim ()];

   for (s = 0; s < config.J; s++) {
      for (int i = 0; i <= config.compon[s].k; i++)
         coef[s][i] = 0;
   }
}


//===========================================================================

void Finalize ()
{
   int s;
   for (s = 0; s < config.J; s++)
      delete[]reg[s];
   delete[]reg;

   for (s = 0; s < config.J; s++)
      delete[]zone[s];
   delete[]zone;

   for (s = 0; s < config.J; s++) {
      if (config.compon[s].PerMax)
         delete primJ[s];
      delete compJ[s];
   }

   delete[]CoNoElemPrim;
   delete[]CoElemPrim;
   delete[]TotEP;
   delete[]CoPerMax;
   delete[]CoRegions;
   delete[]Hk;
   delete[]H;
   delete[]minVal;

   delete[]compJ;
   delete[]primJ;
   delete[]pool;
   delete normal;
}


//===========================================================================

void InitZones ()
{
   for (int s = 0; s < config.J; s++) {
      for (int i = 1; i <= config.compon[s].k; i++) {
         zone[s][i].init (config.compon[s], s, i);
      }
   }
}


//===========================================================================

void SeekGen (int j)
{
   if (j >= config.J) {
      Test ();
      CompareMerit ();
      ++numATried;
      return ;
   }

   for (int region = 0; region < config.compon[j].numReg; region++) {
      if (timer.timeOver (config.duration))
         return ;

      if (config.compon[j].searchMethod == EXHAUST) {
         ExamAllZones (j, config.compon[j].k);
      } else {
         ChoisirBornes (j);
         ExamRegion (j, config.compon[j].k);
      }
   }
}


//===========================================================================

}                                 // namespace


//===========================================================================

int main (int argc, char **argv)
{
   if (readConfigFile (argc, argv))
      return -1;
   config.write ();
   Init ();
   timer.init ();
   Zone::initFrontieres (config);
   InitZones ();
   SeekGen (0);
   SortBestGen ();
   PrintResults ();
   Finalize ();
}
