#ifndef MERIT_H
#define MERIT_H
#include <vector>
#include <string>
#include "latcommon/Const.h"


namespace LatMRG {

/**
 * This class is used to keep the values of test results for a lattice. The
 * values of the figure of merit in each dimension for a given lattice are
 * kept in vectors. For reasons of efficiency, the values kept in merit
 * objects may not always be the merit itself but instead a simple function
 * of the merit. Only at the end of the tests will the real merit be printed.
 * For the spectral test, in the case of the \f${\mathcal{L}}_2\f$ norm, the
 * square of the length of the shortest non-zero vector in the lattice is
 * kept in the merit object. However, in the case of the
 * \f${\mathcal{L}}_1\f$ norm, the length of the shortest non-zero vector in
 * the lattice is kept in the merit object.
 *
 * As an example, the real normalized merit for the spectral test (the
 * shortest vector normalized with the constants for laminated lattices in
 * the case of the \f${\mathcal{L}}_2\f$ norm) for the LCG with \f$m=
 * 2^{31} -1\f$ and \f$a=16807\f$, in dimensions up to 10, is given in the
 * following table:
 *
 * <center>
 *
 * <table class="LatSoft-table LatSoft-has-hlines">
 * <tr class="bt">
 *   <td class="c bl br"></td>
 *   <td colspan="2">\f${\mathcal{L}}_2\f$ norm</td>
 *   <td colspan="2">\f${\mathcal{L}}_1\f$ norm</td>
 * </tr><tr class="bt">
 *   <td class="c bl br">dim</td>
 *   <td class="l bl br">Dual</td>
 *   <td class="l bl br">Primal</td>
 *   <td class="l bl br">Dual</td>
 *   <td class="l bl br">Primal</td>
 * </tr><tr class="bt">
 *   <td class="c bl br">2</td>
 *   <td class="l bl br">0.33751</td>
 *   <td class="l bl br">0.33751</td>
 *   <td class="l bl br">0.25647</td>
 *   <td class="l bl br">0.25647</td>
 * </tr><tr>
 *   <td class="c bl br">3</td>
 *   <td class="l bl br">0.44118</td>
 *   <td class="l bl br">0.54043</td>
 *   <td class="l bl br">0.32637</td>
 *   <td class="l bl br">0.51556</td>
 * </tr><tr>
 *   <td class="c bl br">4</td>
 *   <td class="l bl br">0.57519</td>
 *   <td class="l bl br">0.61619</td>
 *   <td class="l bl br">0.57143</td>
 *   <td class="l bl br">0.55061</td>
 * </tr><tr>
 *   <td class="c bl br">5</td>
 *   <td class="l bl br">0.73612</td>
 *   <td class="l bl br">0.61872</td>
 *   <td class="l bl br">0.67539</td>
 *   <td class="l bl br">0.58963</td>
 * </tr><tr>
 *   <td class="c bl br">6</td>
 *   <td class="l bl br">0.64541</td>
 *   <td class="l bl br">0.5889</td>
 *   <td class="l bl br">0.58879</td>
 *   <td class="l bl br">0.52362</td>
 * </tr><tr>
 *   <td class="c bl br">7</td>
 *   <td class="l bl br">0.57112</td>
 *   <td class="l bl br">0.60365</td>
 *   <td class="l bl br">0.5</td>
 *   <td class="l bl br">0.45958</td>
 * </tr><tr>
 *   <td class="c bl br">8</td>
 *   <td class="l bl br">0.60961</td>
 *   <td class="l bl br">0.44125</td>
 *   <td class="l bl br">0.50909</td>
 *   <td class="l bl br">0.35768</td>
 * </tr><tr>
 *   <td class="c bl br">9</td>
 *   <td class="l bl br">0.57732</td>
 *   <td class="l bl br">0.65822</td>
 *   <td class="l bl br">0.46667</td>
 *   <td class="l bl br">0.58086</td>
 * </tr><tr>
 *   <td class="c bl br">10</td>
 *   <td class="l bl br">0.65033</td>
 *   <td class="l bl br">0.53741</td>
 *   <td class="l bl br">0.5</td>
 *   <td class="l bl br">0.48685</td>
 * </tr>
 * </table>
 *
 * </center>
 *
 */
class Merit {
public:

/**
 * Constructor. This object may contain figures of merit for dimensions up to
 * `maxDim`.
 */
Merit (int maxDim);

   /**
    * Copy Constructor.
    */
   Merit (const Merit & mer);

   /**
    * Assignment operator.
    */
   Merit & operator= (const Merit & mer);

   /**
    * Destructor.
    */
   ~Merit ();

   /**
    * Returns the effective dimension of the vectors.
    */
   int getDim () const { return m_normVal.size() - 1; }

   /**
    * Returns the *unnormalized* value of the merit in dimension \f$j\f$.
    */
   double & getMerit (int j)  {
      // return m_val.at(j);
      return m_val[j];
     }

   /**
    * Returns the *normalized* value of the merit in dimension \f$j\f$.
    */
   double getNormVal (int j) const  {
      // return m_normVal.at(j);
      return m_normVal[j];
     }

   /**
    * Returns the *normalized* value of the merit in dimension \f$j\f$.
    * (Same as `getNormVal`.)
    */
   double & operator[] (int j)  {
      // return m_normVal.at(j);
      return m_normVal[j];
     }

   /**
    * Sets the value of the merit to \f$x\f$ for all dimensions.
    */
   void set (double x);

   /**
    * Returns the smallest *normalized* value of the merit in the
    * dimension range `from` to `T` (inclusive). As a side effect, the
    * dimension where the smallest value occurs is returned in `dimWorst`.
    */
   double getST (int from, int T, int & dimWorst);

   /**
    * Returns the smallest *normalized* value of the merit in the
    * dimension range `from` to `T` (inclusive). As a side effect, the
    * dimension where the smallest value occurs is set in the private
    * variable <tt>m_dimWorst</tt>.
    */
   double getST (int from, int T){return getST(from, T, m_dimWorst);}

   /**
    * Returns the dimension of the worst (the smallest usually) figure of
    * merit.
    */
   int getDimWorst () const { return m_dimWorst; }

   /**
    * Sets the worst figure of merit to \f$x\f$.
    */
   void setWorstMerit (double x) { m_worstMerit = x; }

   /**
    * Returns the worst figure of merit. It must have been calculated
    * before.
    */
   double getWorstMerit () const { return m_worstMerit; }

   /**
    * Returns both the *unnormalized* and the *normalized* values of the
    * merit from dimension `from` to `T` as a string. If `rac` is `true`,
    * the square root of the values of the merit is returned; otherwise,
    * the values themselves are returned. If `invert` is `true`, the
    * inverses of the length of the shortest vector are returned;
    * otherwise, the lengths are returned.
    */
   std::string toString (int from, int T, bool rac, bool invert) const;
private:

/**
 * The *unnormalized* values of the figure of merit for each dimension.
 */
std::vector<double> m_val;

   /**
    * The *normalized* values of the figure of merit. The corresponding
    * values of <tt>m_val</tt> are normalized so that <tt>m_normVal</tt>
    * usually takes values in [0, 1].
    */
   std::vector<double> m_normVal;

   /**
    * The worst value of the normalized figure of merit for this object.
    */
   double m_worstMerit;

   /**
    * The dimension where the worst value of the merit occurs.
    */
   int m_dimWorst;
};

}
#endif
