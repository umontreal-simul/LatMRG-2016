#ifndef SEEK_CONFIG_H
#define SEEK_CONFIG_H
#include <string>
#include <vector>

#include "latcommon/Types.h"
#include "latcommon/Const.h"
#include "Modulus.h"
#include "latcommon/Util.h"


namespace LatMRG {

struct Component {
   LatCommon::GenType genType;        // Generator type: MRG, MWC
   Modulus modulus;        // Modulus m
   int k;                  // Generator order
   bool PerMax;            // True if maximal period is required, else false
   LatCommon::ImplemCond implemCond;
   int NumBits;
   int HighestBit;
   int ncoef;
   int *Icoef;
   LatCommon::DecompType F1; 
   std::string file1;
   LatCommon::DecompType F2; 
   std::string file2;
   LatCommon::SearchMethod searchMethod; 
   int numReg, H, Hk;
   MVect b, c;             // intervals where to search the coefs 
   bool ApproxTotGen;
};

/**
 * \f$\clubsuit\f$ Compléter la description de la classe
 *
 * <strong> À EXAMINER:</strong> il pourrrait y avoir des avantages à inclure
 * une variable `LatConfig` à l’intérieur de `SeekConfig`. Plusieurs
 * variables n’auraient pas à être répétées dans `SeekConfig`, et la
 * recherche appelle les tests avec des paramètres de `LatConfig`, parfois
 * explicitement, ce qui cause des problèmes (cas PALPHA). Cela pourrait
 * simplifier le code.
 *
 */
class SeekConfig {
public:
   static const int MAX_WORD_SIZE = 64;
   SeekConfig();
   ~SeekConfig();
   void write();
   MScal* getMs() { return Ms; }
   int* getKs()  { return Ks; }
   int getMaxK();
   int getMaxTd();       
   bool readGenFile; // read generator from file
   std::string fileName;  // file name
   int J;            // nombre de generateurs dans la combinason
   std::vector<Component> compon;
   MScal* Ms;
   MVect* As;
   int* Ks;
   int C; 
   double* minMerit; 
   double* maxMerit; // C values
   int* numGen; // C values
LatCommon::CriterionType criter; 
   LatCommon::NormaType normaType;

/**
 * The number of series of projections (see `td` below). The classical case
 * corresponds to \f$d=1\f$, for which the chosen test is done on all
 * successive dimensions from `td[0] = fromDim`, up to and including `td[1] =
 * toDim`.
 */
// Criterion <Norm> 

   int d;

   /**
    * Array containing the maximal dimensions for the projections in each
    * category. `td[1]` is the maximal dimension for successive
    * 1-dimensional projections, `td[2]` is the maximal dimension for
    * 2-dimensional projections, `td[3]` is the maximal dimension for
    * 3-dimensional projections, and so on. However, the value of `td[0]`
    * is the minimal dimension for the case of successive dimensions.
    */
   int *td;
int getMaxDim() { return td[1]; }
   LatCommon::LatticeType latType;

/**
 * If this flag is `true`, the test is to be applied on the dual lattice. If
 * it is `false`, the test is to be applied on the primal lattice.
 */
bool dualF;

   /**
    * If `invertF` is `true`, the inverse of the length of the shortest
    * vector will be printed in the results. Otherwise, the length itself
    * will be printed.
    */
   bool invertF;

   /**
    * The value of \f$\alpha\f$ for the \f$P_{\alpha}\f$ test.
    */
   int alpha;
   int lacGroupSize;
   int lacSpacing;
   long maxNodesBB;
   double duration;
   long seed;  // seed of the random number generator
   LatCommon::OutputType outputType;
};

}
#endif // SEEK_CONFIG_H
