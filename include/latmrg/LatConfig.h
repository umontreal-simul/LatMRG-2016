#ifndef LATCONFIG_H
#define LATCONFIG_H
#include "latcommon/Types.h"
#include "latcommon/Const.h"
#include "latmrg/MRGComponent.h"


namespace LatMRG {

/**
 * This class is used to save the configuration of a lattice test. It is used
 * to keep all the parameters read in the data file and passed to different
 * methods for the spectral, Beyer or \f$P_{\alpha}\f$ tests.
 *
 */
class LatConfig {
public:
   static const int MAX_WORD_SIZE = 64;

   /**
    * Constructor.
    */
   LatConfig();

   /**
    * Destructor.
    */
   ~LatConfig();

   /**
    * Frees the memory used by this object.
    */
   void kill();

   /**
    * For debugging: writes the configuration on the console.
    */
   void write();

   /**
    * Reinitializes this object and allocates enough memory for \f$j\f$
    * MRGs.
    */
   void setJ (int j);

   /**
    * This flag is set `true` if the generator is to be read from a file,
    * otherwise it is set `false`.
    */
   bool readGenFile;

   /**
    * If `readGenFile` is `true`, the name of the file from which to read
    * the generator.
    */
   std::string fileName;

   /**
    * The number of MRG components.
    */
   int J;

   /**
    * The array of generator types for each MRG component. See module
    * `Const` for more details.
    */
   LatCommon::GenType *genType;

   /**
    * The array of MRG components which describe the combined generator.
    */
   MRGComponent **comp;

   /*
    * The minimal dimension for which the test will be performed. ************
    * REMPLACÉ PAR `td[0]`. À ÉLIMINER.
    */
   //   int fromDim;

   /*
    * The maximal dimension for which the test will be performed. ************
    * REMPLACÉ PAR `td[1]`. À ÉLIMINER.
    */
   //   int toDim;

   /**
    * The number of categories of projections (see `td` below). The
    * classical case corresponds to \f$d=1\f$, for which the chosen test
    * is done on all successive dimensions from `td[0] = fromDim`, up to
    * and including `td[1] = toDim`.
    */
   int d;

   /**
    * Array containing the maximal dimensions for the projections in each
    * category. `td[1]` is the maximal dimension for successive
    * 1-dimensional projections, `td[2]` is the maximal dimension for
    * 2-dimensional projections, `td[3]` is the maximal dimension for
    * 3-dimensional projections, and so on. The value of `td[0]` is the
    * minimal dimension for the case of successive dimensions.
    */
   int *td;

   /**
    * The criterion for which the test will be performed. See module
    * `Const` for the possible criterion types.
    */
   LatCommon::CriterionType criter;

   /**
    * The bound used for the normalization in the definition of \f$S_t\f$.
    * Only applicable for the spectral test.
    */
   LatCommon::NormaType norma;

   /**
    * The type of the calculation in the \f$P_{\alpha}\f$ test.
    */
   LatCommon::CalcType calcPalpha;

   /**
    * The value of \f$\alpha\f$ for the \f$P_{\alpha}\f$ test.
    */
   int alpha;

   /**
    * The values of \f$B_i\f$, \f$i=0, 1,…,\f$<tt>toDim</tt> for the
    * \f$P_{\alpha}\f$ test.
    */
   std::vector<double> Beta;

   /**
    * The seed for the generator in the \f$P_{\alpha}\f$ test.
    */
   int seed;

   /**
    * This flag is set `true` if the test is to be applied on the dual
    * lattice. If it is `false`, the test is applied on the primal
    * lattice.
    */
   bool dualF;

   /**
    * If `invertF` is `true`, the inverse of the length of the shortest
    * vector will be printed in the results. Otherwise, the length itself
    * will be printed.
    */
   bool invertF;

   /**
    * This flag indicates to print more detailed results if `detailF`
    * \f$>0\f$. Default value: 0.
    */
   int detailF;

   /**
    * Norm used to measure the length of vectors. See module `Const` for a
    * definition of the possible norms.
    */
   LatCommon::NormType norm;

   /**
    * Indicates the type of lattice used in the test. See `Const` for a
    * definition of the possible lattice types.
    */
   LatCommon::LatticeType latType;

   /**
    * This flag is set `true` if the test is applied for lacunary indices.
    * If it is `false`, the test is applied for successive indices.
    */
   bool lacunary;

   /**
    * Used for lacunary indices. If the respective values are \f$s\f$ and
    * \f$d\f$, then the program will analyze the lattice structure of vectors
    * formed by groups of \f$s\f$ successive values, taken \f$d\f$ values apart,
    * i.e. groups of the form \f$(u_{i+1}, …, u_{i+s}, u_{i+d+1}, …, u_{i+d+s},
    * u_{i+2d+1}, …, u_{i+2d+s}, …)\f$.
    */
   int lacGroupSize;

   /**
    * \copydoc lacGroupSize
    */
   NTL::ZZ lacSpacing;

   /**
    * The lacunary indices, either read explicitly or computed from
    * `lacGroupSize` and `lacSpacing`.
    */
   BVect Lac;

   /**
    * Is `true` when the modulus of congruence \f$m\f$ is a prime number,
    * is `false` otherwise.
    */
   bool primeM;

   /**
    * If `true`, the program will verify that the modulus \f$m\f$ is a
    * prime number. If `false`, will not verify it.
    */
   bool verifyM;

   /**
    * Is `true` when the generator has maximal period, is `false`
    * otherwise.
    */
   bool maxPeriod;

   /**
    * If `true`, the program will verify that the generator has maximal
    * period. If `false`, will not verify it.
    */
   bool verifyP;

   /**
    * The maximum number of nodes to be examined in any given
    * branch-and-bound procedure when computing \f$d_t\f$ or \f$q_t\f$.
    */
   long maxNodesBB;

   /**
    * File format used to store the results. See `Const` for a definition
    * of the possible output types.
    */
   LatCommon::OutputType outputType;
};

}
#endif
