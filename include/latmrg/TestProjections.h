#ifndef TESTPROJECTIONS_H
#define TESTPROJECTIONS_H
#include "latmrg/IntLattice.h"
#include "LatticeTest.h"
#include "latcommon/Weights.h"
#include "ProjIterator.h"
#include "latcommon/CoordinateSets.h"
#include "Writer.h"


namespace LatMRG {

/**
 * Implements methods used to calculate the worst-case figure of merit,
 * defined as follows:
 * \f[
 *   M_{t_1,…,t_d} = \min\left[ \min_{k+1\le t\le t_1} \frac{\ell_t}{\ell_t^*(m^k)},\; \min_{2\le s\le k}\; \min_{I\in S(s,t_s)} \frac{\ell_I}{m}, \; \min_{k+1\le s\le d} \;\min_{I\in S(s,t_s)} \frac{\ell_I}{\ell_s^*(m^k)} \right],
 * \f]
 * where \f$S(s, t_s) = \left\{I = \{i_1,…,i_s\} \mid1 = i_1 < \cdots< i_s
 * \le t_s\right\}\f$. In other words, this figure of merit applies the
 * chosen test over the \f$s\f$ successive dimensions for all \f$s
 * \le t_1\f$, and over the \f$r\f$ nonsuccessive dimensions of the lattice
 * of dimension \f$t_r\f$ for all \f$2 \le r \le d\f$. For *dimension
 * stationary* lattices, for example Korobov lattices, only the sets of
 * dimensions whose first coordinate is 1 need to be considered.
 *
 * Here is an example of the spectral test with different projections for a
 * Korobov lattice with \f$m=1021\f$ and \f$a=333\f$, when we consider the
 * dual lattice with normalization BESTLAT. The figure of merit is
 * \f$M_{10,8,6,5}\f$ and the resulting minimal merit is 0.440728, for
 * projection 1,6. The length is the length of the shortest vector of the
 * lattice with the \f${\mathcal{L}}_2\f$ norm.
 *
 * <center>
 *
 * <table class="LatSoft-table LatSoft-has-hlines">
 * <tr class="bt">
 *   <td class="l">projections</td>
 *   <td class="l">length</td>
 *   <td class="l">merit</td>
 * </tr><tr class="bt">
 *   <td class="l">2</td>
 *   <td class="l">22.2036</td>
 *   <td class="l">0.64666</td>
 * </tr><tr>
 *   <td class="l">3</td>
 *   <td class="l">7</td>
 *   <td class="l">0.619324</td>
 * </tr><tr>
 *   <td class="l">4</td>
 *   <td class="l">3.87298</td>
 *   <td class="l">0.576145</td>
 * </tr><tr>
 *   <td class="l">5</td>
 *   <td class="l">3.74166</td>
 *   <td class="l">0.760239</td>
 * </tr><tr>
 *   <td class="l">6</td>
 *   <td class="l">2</td>
 *   <td class="l">0.488395</td>
 * </tr><tr>
 *   <td class="l">7</td>
 *   <td class="l">2</td>
 *   <td class="l">0.552276</td>
 * </tr><tr>
 *   <td class="l">8</td>
 *   <td class="l">2</td>
 *   <td class="l">0.594822</td>
 * </tr><tr>
 *   <td class="l">9</td>
 *   <td class="l">2</td>
 *   <td class="l">0.654906</td>
 * </tr><tr>
 *   <td class="l">10</td>
 *   <td class="l">2</td>
 *   <td class="l">0.697213</td>
 * </tr><tr class="bt">
 *   <td class="l">1,3</td>
 *   <td class="l">25.4951</td>
 *   <td class="l">0.742522</td>
 * </tr><tr>
 *   <td class="l">1,4</td>
 *   <td class="l">20.6155</td>
 *   <td class="l">0.600409</td>
 * </tr><tr>
 *   <td class="l">1,5</td>
 *   <td class="l">30.6105</td>
 *   <td class="l">0.891502</td>
 * </tr><tr>
 *   <td class="l">1,6</td>
 *   <td class="l">15.1327</td>
 *   <td class="l">0.440728</td>
 * </tr><tr>
 *   <td class="l">1,7</td>
 *   <td class="l">16.6433</td>
 *   <td class="l">0.484722</td>
 * </tr><tr>
 *   <td class="l">1,8</td>
 *   <td class="l">32.3883</td>
 *   <td class="l">0.943279</td>
 * </tr><tr class="bt">
 *   <td class="l">1,2,4</td>
 *   <td class="l">7.07107</td>
 *   <td class="l">0.625612</td>
 * </tr><tr>
 *   <td class="l">1,2,5</td>
 *   <td class="l">8.12404</td>
 *   <td class="l">0.718773</td>
 * </tr><tr>
 *   <td class="l">1,2,6</td>
 *   <td class="l">9.48683</td>
 *   <td class="l">0.839346</td>
 * </tr><tr>
 *   <td class="l">1,3,4</td>
 *   <td class="l">9.43398</td>
 *   <td class="l">0.83467</td>
 * </tr><tr>
 *   <td class="l">1,3,5</td>
 *   <td class="l">9.69536</td>
 *   <td class="l">0.857795</td>
 * </tr><tr>
 *   <td class="l">1,3,6</td>
 *   <td class="l">8.12404</td>
 *   <td class="l">0.718773</td>
 * </tr><tr>
 *   <td class="l">1,4,5</td>
 *   <td class="l">7.54983</td>
 *   <td class="l">0.66797</td>
 * </tr><tr>
 *   <td class="l">1,4,6</td>
 *   <td class="l">8.12404</td>
 *   <td class="l">0.718773</td>
 * </tr><tr>
 *   <td class="l">1,5,6</td>
 *   <td class="l">6.78233</td>
 *   <td class="l">0.600066</td>
 * </tr><tr class="bt">
 *   <td class="l">1,2,3,5</td>
 *   <td class="l">5.74456</td>
 *   <td class="l">0.854561</td>
 * </tr><tr>
 *   <td class="l">1,2,4,5</td>
 *   <td class="l">3.87298</td>
 *   <td class="l">0.576145</td>
 * </tr><tr>
 *   <td class="l">1,3,4,5</td>
 *   <td class="l">5.47723</td>
 *   <td class="l">0.814792</td>
 * </tr>
 * </table>
 *
 * </center>
 *
 */
class TestProjections {
public:

   /**
    * `master` is the original lattice for which we want to calculate the
    * \f$M\f$ merit as defined above. `lattice` is the working lattice used for
    * intermediate calculations; all the projections from `master` are stored in
    * `lattice` before calling the test. `test` is the lattice test to be
    * applied on the different projections. `d` is the last element of array
    * `td`, which gives the maximal dimensions for the different projections.
    * `td[0]` and `td[1]` gives the minimal and the maximal dimensions for the
    * test applied on successive dimensions. For example, if `d` = 3 and `td` =
    * [2, 32, 16, 12], then the test will be applied on all successive
    * dimensions \f$2 \le t \le32\f$, on all possible 2-dimensional projections
    * for dimensions \f$\le16\f$, and on all possible 3-dimensional projections
    * for dimensions \f$\le12\f$. The results of the tests will be outputted on
    * `Writer`.
    */
   TestProjections (IntLattice *master, IntLattice *lattice,
                       LatticeTest *test, int td[], int d);

   /**
    * Destructor.
    */
   ~TestProjections ();

   /**
    * Sends the output to `rw`. If this method is not called, output will
    * be sent to standard output.
    */
   void setOutput (Writer * rw);

   /**
    * If `flag` is `true`, the tests will be applied on the *dual*
    * lattice; if it is `false`, the tests will be applied on the *primal*
    * lattice.
    */
   void setDualFlag (bool flag);

   /**
    * Sets the value of the <tt>m_invertF</tt> flag. If `invertF` is
    * `true`, the inverse of the length of the shortest vector will be
    * printed in the results. Otherwise, the length itself will be
    * printed.
    */
   void setInvertFlag (bool flag);

   /**
    * If flag is `true`, the value of the merit will be printed for all
    * projections, otherwise not. This flag should be set `false` in the
    * case of the `seek*` programs, and `true` otherwise.
    */
   void setPrintF (bool flag);

   /**
    * Builds the basis (and dual basis) of the projection `proj` for this
    * lattice. The result is placed in the work lattice. The basis is
    * triangularized to form a proper basis. For example, if \f$d=2\f$ and
    * `indices` = \f$[1, 4]\f$, then a 2-dimensional basis is built using
    * coordinates 1 and 4 of the master basis.
    */
   void build (const LatCommon::Coordinates & proj);

   /**
    * Calculates the \f$M_{t_1, …, t_d}\f$ merit by running all the tests
    * for the different projections. If set `true`, `stationary` means
    * that the lattice is known to be *dimension stationary*. If set
    * `true`, the flag `last` means that only the projections which
    * include the last dimension of the lattice will be considered. This
    * is good for performing incremental searches for good lattices.
    * `minVal` is the minimal value of the normalized merit in each
    * dimension for a lattice to be considered. This method returns the
    * worst value of the merit over all projections.
    */
   double run (bool stationary, bool last, double minVal[]);

   /**
    * As method `run` above, but with the weights `weights`.
    */
   double run (bool stationary, bool last, double minVal[],
               const LatCommon::Weights & weights);

   /**
    * Calculates the number of projections, given the parameters of this
    * object. The flag `stationary` must be set `true` if the lattice is
    * *dimension stationary*. If set `true`, the flag `last` means that
    * only the projections which include the last dimension of the lattice
    * will be considered.
    */
   int calcNumProjections (bool stationary, bool last);

   /**
    * Returns the number of projections considered for the last call to
    * `run`. If the test was interrupted because the lattice was rejected
    * as uninteresting, the number returned will be less than the total
    * number of projections.
    */
   int getNumProjections () { return m_numproj; }

protected:

   /**
    * Run only for projections given by the iterator `projit`.
    */
   double run (ProjIterator & projit, double minVal[],
                  const LatCommon::Weights & weights);

   /**
    * Lattice on which the \f$M_{t_1, …, t_d}\f$ merit will be calculated.
    */
   IntLattice* m_master;

   /**
    * Working lattice which is used to test the different projections.
    */
   IntLattice* m_lattice;

   /**
    * The lattice test used to calculate the merit.
    */
   LatticeTest* m_test;

   /**
    * Weights.
    */
   double* m_weightsTemp;

   /**
    * If `true`, the test is applied on the dual lattice, otherwise on the
    * primal lattice.
    */
   bool m_dualF;

   /**
    * If `true`, the inverse length of the shortest vector is printed for
    * all projections, otherwise the length is printed.
    */
   bool m_invertF;

   /**
    * If `true`, the value of the merits is printed for all projections,
    * otherwise not.
    */
   bool m_printF;

   /**
    * If `true`, the square root of the values of the merit is printed
    * after the test; otherwise the values themselves are printed.
    */
   bool m_racF;

   /**
    * The maximal dimensions for each kind of projections.
    * <tt>m_td[0]</tt> is the minimal dimension for successive dimensions.
    * <tt>m_td[1]</tt> is the maximal dimension for successive dimensions
    * (1-dimensional projections), <tt>m_td[2]</tt> is the maximal
    * dimension for 2-dimensional projections, <tt>m_td[3]</tt> is the
    * maximal dimension for 3-dimensional projections, and so on. The last
    * element, <tt>m_td[m_d]</tt>, is the maximal dimension for
    * <tt>m_d</tt>-dimensional projections.
    */
   int *m_td;

   /**
    * The number of kinds of projections. Also the number of elements of
    * array <tt>m_td</tt> is <tt>m_d + 1</tt>.
    */
   int m_d;

   /**
    * The number of projections.
    */
   int m_numproj;
private:

/**
 * Output will be written on <tt>m_writer</tt>. By default, it is written on
 * standard output.
 */
Writer *m_writer;
   bool m_wrFlag;
};

}
#endif
