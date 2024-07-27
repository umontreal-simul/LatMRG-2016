#ifndef PARAMREADER_H
#define PARAMREADER_H
#include "NTL/ZZ.h"
#include "latcommon/Types.h"
#include "latcommon/Util.h"
#include "latcommon/Const.h"
#include "latmrg/MRGComponent.h"
#include <string>
#include <vector>


namespace LatMRG {

/**
 * Utility class used to read basic parameter fields in a configuration file.
 * Lines whose first non-blank character is a <tt>#</tt> are considered as
 * comments and discarded.
 *
 */
class ParamReader {
public:
   static const int MAX_WORD_SIZE = 64;

/**
 * Constructor.
 */
ParamReader();

   /**
    * Constructor. Opens the file `fileName`.
    */
   ParamReader (std::string fileName);

   /**
    * Destructor.
    */
   ~ParamReader();

   /**
    * Reads all the lines from the file and stores them into this object’s
    * buffer. Lines whose first non-blank character is a <tt>#</tt> are
    * considered as comments and discarded. Empty lines are also
    * discarded.
    */
   void getLines();

   /**
    * Puts into `field` the <tt>pos</tt>-th string token from line `ln`.
    */
   void getToken (std::string & field, unsigned int ln, unsigned int pos);

   /**
    * Splits line `ln` from the file into several string tokens. Separator
    * characters are defined in function `IsDelim`. Tokens are stored in
    * vector `tokens`.
    */
   int tokenize (std::vector<std::string> & tokens, unsigned int ln);

   /**
    * Reads a string from the <tt>pos</tt>-th token of the <tt>ln</tt>-th
    * line into `field`.
    */
   void readString (std::string & field, unsigned int ln, unsigned int pos);

   /**
    * Reads a boolean from the <tt>pos</tt>-th token of the <tt>ln</tt>-th
    * line into `field`.
    */
   void readBool (bool & field, unsigned int ln, unsigned int pos);

   /**
    * Reads a character from the <tt>pos</tt>-th token of the
    * <tt>ln</tt>-th line into `field`.
    */
   void readChar (char & field, unsigned int ln, unsigned int pos);

   /**
    * Reads an integer from the <tt>pos</tt>-th token of the
    * <tt>ln</tt>-th line into `field`.
    */
   void readInt (int & field, unsigned int ln, unsigned int pos);

   /**
    * Reads a long from the <tt>pos</tt>-th token of the <tt>ln</tt>-th
    * line into `field`.
    */
   void readLong (long & field, unsigned int ln, unsigned int pos);

   /**
    * Reads a large integer from the <tt>pos</tt>-th token of the
    * <tt>ln</tt>-th line into `field`.
    */
   void readZZ (NTL::ZZ & field, unsigned int ln, int pos);

   /**
    * Reads a double from the <tt>pos</tt>-th token of the <tt>ln</tt>-th
    * line into `field`.
    */
   void readDouble (double & field, unsigned int ln, unsigned int pos);

   /**
    * Reads a `GenType` from the <tt>pos</tt>-th token of the
    * <tt>ln</tt>-th line into `field`.
    */
   void readGenType (LatCommon::GenType & field, unsigned int ln,
                     unsigned int pos);

   /**
    * Reads \f$b\f$, \f$e\f$ and \f$c\f$, starting at the <tt>pos</tt>-th
    * token of the <tt>ln</tt>-th line and uses them to define \f$r\f$.
    * The numbers in the data file may be given in one of the two
    * following formats:
    *
    * \f$\bullet\f$ A single integer giving the value of \f$r=b\f$
    * directly on a line. In that case, one sets \f$e=c=0\f$. <br>
    * \f$\bullet\f$ Three integers \f$b\f$, \f$e\f$, \f$c\f$ on the same
    * line, separated by at least one blank. The \f$r\f$ value will be set
    * as \f$r=b^e+c\f$ if \f$b>0\f$, and \f$r= -(|b|^e+c)\f$ if \f$b<0\f$.
    * One must have \f$e\ge0\f$. For example, \f$(b, e, c) = (2,
    * 5, -1)\f$ will give \f$r=31\f$, while \f$(b, e, c) = (-2, 5, -1)\f$
    * will give \f$r=-31\f$.
    */
   void readNumber3 (MScal & r, long & b, long & e, long & c,
                     unsigned int ln, unsigned int pos);

   /**
    * Reads a `BScal` from the <tt>pos</tt>-th token of the <tt>ln</tt>-th
    * line into `field`.
    */
   void readBScal (BScal & field, unsigned int ln, int pos);

   /**
    * Reads a `MScal` from the <tt>pos</tt>-th token of the <tt>ln</tt>-th
    * line into `field`.
    */
   void readMScal (MScal & field, unsigned int ln, unsigned int pos);

   /**
    * Reads a `RScal` from the <tt>pos</tt>-th token of the <tt>ln</tt>-th
    * line into `field`.
    */
   void readRScal (RScal & field, unsigned int ln, unsigned int pos);

   /**
    * Reads `num` tokens (from the <tt>pos</tt>-th token of the
    * <tt>ln</tt>-th line) into `field`, starting at index \f$j\f$ of
    * `field`.
    */
   void readMVect (MVect & field, unsigned int & ln, unsigned int pos,
                   unsigned int num, int j);

   /**
    * Reads `num` tokens (from the <tt>pos</tt>-th token of the
    * <tt>ln</tt>-th line) into `field`, starting at index \f$j\f$ of
    * `field`.
    */
   void readIntVect (int* field, unsigned int ln, unsigned int pos,
                     unsigned int num, int j);

   /**
    * Reads `num` tokens (from the <tt>pos</tt>-th token of the
    * <tt>ln</tt>-th line) into `field`, starting at index \f$j\f$ of
    * `field`.
    */
   void readDoubleVect (double* field, unsigned int ln, unsigned int pos,
                        unsigned int num, int j);

   /**
    * Reads `2k MScal` tokens into vectors `B` and `C`, starting at the
    * <tt>ln</tt>-th line. These represent a box \f$[B_i, C_i]\f$, \f$i =
    * 1, 2, …, k\f$. The \f$B_i, C_i\f$ must be given in the order \f$B_1,
    * C_1, B_2, C_2, …, B_k, C_k\f$, each on a line of its own. Each
    * coefficient may be given in the form described in `readNumber3`
    * above.
    */
   void readInterval (MVect & B, MVect & C, unsigned int & ln, int k);

   /**
    * Reads a criterion from the <tt>pos</tt>-th token of the
    * <tt>ln</tt>-th line into `field`.
    */
   void readCriterionType (LatCommon::CriterionType & field, unsigned int ln,
                           unsigned int pos);

   /**
    * Reads a norm from the <tt>pos</tt>-th token of the <tt>ln</tt>-th
    * line into `field`.
    */
   void readNormType (LatCommon::NormType & field, unsigned int ln,
                      unsigned int pos);

   /**
    * Reads a type of normalization from the <tt>pos</tt>-th token of the
    * <tt>ln</tt>-th line into `field`.
    */
   void readNormaType (LatCommon::NormaType & field, unsigned int ln,
                       unsigned int pos);

   /**
    * Reads the type of calculation to do (<tt>PAL</tt> or <tt>BAL</tt>)
    * for the `PALPHA` test from the <tt>pos</tt>-th token of the
    * <tt>ln</tt>-th line into `field`.
    */
   void readCalcType (LatCommon::CalcType & field, unsigned int line,
                      unsigned int pos);

   /**
    * Reads the decomposition type from the <tt>pos</tt>-th token of the
    * <tt>ln</tt>-th line into `field`.
    */
   void readDecompType (LatCommon::DecompType & field, unsigned int line,
                        unsigned int pos);

   /**
    * Reads a lattice type from the <tt>pos</tt>-th token of the
    * <tt>ln</tt>-th line into `field`.
    */
   void readLatticeType (LatCommon::LatticeType & field, unsigned int ln,
                         unsigned int pos);

   /**
    * Reads the initial state for an orbit of a combined generator with
    * \f$J\f$ components `comp`, starting at the <tt>ln</tt>-th line. Each
    * line must have a set of \f$k_j\f$ numbers which are the seeds of the
    * orbit for each component, where \f$k_j\f$ is the order of the
    * \f$j\f$-th component. On exiting this method, the line number is
    * reset to `ln` + \f$J\f$.
    */
   void readOrbit (int J, MRGComponent **comp, unsigned int & ln);

   /**
    * In the case where `lat = PRIMEPOWER`, checks that the modulus
    * \f$m\f$ is given in the form \f$m=b^e\f$, (that is \f$e>0\f$ and
    * \f$c=0\f$. Checks also that the order \f$k=1\f$. If these conditions
    * are not satisfied, stops the program.
    */
   bool checkPrimePower (LatCommon::LatticeType lat, long e, long c, int k);

   /**
    * Reads the fields, starting at the <tt>ln</tt>-th line, for a
    * generator of order \f$k\f$ and for a test in dimensions `fromDim` to
    * `toDim`, to determine the lacunary indices, which are positive
    * integers that will be read into `Lac`. The vector of lacunary
    * indices `Lac` is created here with the appropriate dimension. If the
    * first two values read are \f$s = \f$ `lacGroupSize` and \f$d = \f$
    * `lacSpacing` and if \f$s>0\f$, then what will be analyzed is the
    * lattice structure of vectors of the form \f$(u_{i+1}, …, u_{i+s},
    * u_{i+d+1},…, u_{i+d+s}, u_{i+2d+1},…, u_{i+2d+s}, …)\f$, formed by
    * groups of \f$s\f$ successive values, taken \f$d\f$ values apart. To
    * analyze lacunary indices that are not evenly spaced, put \f$s
    * = -t\f$ where <em>\f$t = \f$ MaxDim</em> and, on the \f$t\f$ lines
    * that follow, give the \f$t\f$ lacunary indices \f$i_1,…,i_t\f$,
    * which are to be interpreted as in Section  {@link
    * REF__sec1_sec_lacunary lacunary}. In all these cases, the `lacunary`
    * flag is set `true`. To analyze vectors of successive values (the
    * non-lacunary case), take \f$s=d=1\f$ or \f$s\ge\f$ *MaxDim*. In
    * this case, the `lacunary` flag is set `false`.
    */
   void readLacunary (int k, int fromDim, int toDim, unsigned int & ln,
                      bool & lacunary, int & lacGroupSize, NTL::ZZ & lacSpacing,
                      BVect & Lac, LatCommon::GenType genType);

   /**
    * Reads an output form from the <tt>pos</tt>-th token of the
    * <tt>ln</tt>-th line into `field`.
    */
   void readOutputType (LatCommon::OutputType & field, unsigned int ln,
                        unsigned int pos);

   /**
    * Reads an implementation condition from the <tt>pos</tt>-th token of
    * the <tt>ln</tt>-th line into `field`.
    */
   void readImplemCond (LatCommon::ImplemCond & field, unsigned int ln,
                        unsigned int pos);

   /**
    * Reads a search method from the <tt>pos</tt>-th token of the
    * <tt>ln</tt>-th line into `field`.
    */
   void readSearchMethod (LatCommon::SearchMethod & field, unsigned int ln,
                          unsigned int pos);

   /**
    * Checks that the components of \f$A\f$ satisfy \f$-m < A_i < m\f$,
    * for \f$i=1, 2, …, k\f$.
    */
   bool checkBound (const MScal & m, const MVect & A, int k);
private:

/**
 * Internal line buffer.
 */
std::vector<std::string> m_lines;

   /**
    * The path of the opened file.
    */
   std::string m_fileName;

   /**
    * Checks if the character `c` is to be considered as a token separator
    * or not.
    */
   bool isDelim (char c);

   /**
    * Does nothing for now.
    */
   void init() {}
};

}
#endif
