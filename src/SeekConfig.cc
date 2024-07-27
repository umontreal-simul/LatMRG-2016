#include "latmrg/SeekConfig.h"

using namespace std;
using namespace LatCommon;


namespace LatMRG
{

SeekConfig::SeekConfig ()
{
   fileName.reserve(MAX_WORD_SIZE);
   td = new int[1 + d]; //
   C = 1;   // Number of categories; is 1 for now.
   minMerit = new double[C];
   maxMerit = new double[C];
   numGen = new int[C];
}

SeekConfig::~SeekConfig()
{
   fileName.clear();
   for (int i = 0;i < J;i++) {
      compon[i].file1.clear();
      compon[i].file2.clear();
      compon[i].b.kill();
      compon[i].c.kill();
   }
   delete[] td;
   delete[] minMerit;
   delete[] maxMerit;
   delete[] numGen;
}

int SeekConfig::getMaxK()
{
   int max = Ks[0];
   for (int i = 1; i < J; i++) {
      if (Ks[i] > max)
         max = Ks[i];
   }
   return max;
}


int SeekConfig::getMaxTd()
{
   int maxVal = 0;
   for (int i = 0; i < d; i++) {
      maxVal = max(maxVal, td[i]);
   }

   return maxVal;

}

void SeekConfig::write()
{
   cout << "readGenFile: " << boolalpha << readGenFile << endl;
   if (readGenFile)
      cout << "fileName: " << fileName << endl;
   cout << "J: " << J << endl;
   for (int i = 0; i < J; i++) {
      cout << "\n================ Component " << i+1 << " =================\n";
      cout << "typeGen: " << toStringGen (compon[i].genType) << endl;
      cout << "m: " << compon[i].modulus.m << endl;
      cout << "k: " << compon[i].k << endl;
      cout << "PerMax: " << compon[i].PerMax << endl;
      cout << "implemCond: " << toStringImplemCond (compon[i].implemCond) << endl;
      if (POWER_TWO == compon[i].implemCond) {
         cout << "NumBits: " << compon[i].NumBits << endl;
         cout << "HighestBit: " << compon[i].HighestBit << endl;
      }
      cout << "F1: " << toStringDecomp (compon[i].F1) << endl;
      cout << "file1: " << compon[i].file1 << endl;
      cout << "F2: " << toStringDecomp (compon[i].F2) << endl;
      cout << "file2: " << compon[i].file2 << endl;
      cout << "searchMethod: " << toStringSearchMethod (compon[i].searchMethod)
                   << endl;
      if (RANDOM == compon[i].searchMethod) {
         cout << "numReg: " << compon[i].numReg << endl;
         cout << "H: " << compon[i].H << endl;
         cout << "Hk: " << compon[i].Hk << endl;
      }
      cout << "b: " << toString<MVect>(compon[i].b, compon[i].k) << endl;
      cout << "c: " << toString<MVect>(compon[i].c, compon[i].k) << endl;
   }
   cout << "========================= End  ======================\n" << endl;
   cout << "C : " << C << endl;
   cout << "MinMerit: " << toString<double*>(minMerit, 0, C-1) << endl;
   cout << "MaxMerit: " << toString<double*>(maxMerit, 0, C-1) << endl;
   cout << "NumGen: " << toString<int*>(numGen, 0, C-1) << endl;
   cout << "td:   " << toString (td, 0, d) << endl;
   cout << "criter: " << toStringCriterion(criter) << endl;
   cout << "normaType: " << toStringNorma (normaType) << endl;
   if (dualF)
      cout << "lattice: DUAL" << endl;
   else
      cout << "lattice: PRIMAL" << endl;
   cout << "latticeType: " << toStringLattice (latType) << endl;
   cout << "lacGroupSize: " << lacGroupSize << endl;
   cout << "lacSpacing: " << lacSpacing << endl;
   cout << "maxNodesBB: " << maxNodesBB << endl;
   cout << "duration: " << duration << endl;
   cout << "seed: " << seed << endl;
//   cout << "S2: " << s2 << endl;
   cout << "outputType: " << toStringOutput (outputType) << endl;
}

} // namespace
