#ifndef DISCREPANCY_H
#define DISCREPANCY_H
#include "FigureOfMerit.h"
#include "latcommon/ProductWeights.h"
#include "latcommon/OrderDependentWeights.h"
#include <string>


namespace LatMRG {

/**
 * This abstract class is the base class of all discrepancy measures and
 * bounds on discrepancies. The derived classes must implement the pure
 * virtual method `compute`, which computes the discrepancy measure according
 * to its specific definition. The discrepancy is a figure of merit which
 * allows us to measure the "goodness" of a lattice rule in the sense of
 * having a "low discrepancy". This figure of merit may be an upper bound on
 * a specific discrepancy or sometimes, a worst-case error.
 *
 * \remark This class uses functions from MyLib, part of TestU01.
 */
class Discrepancy : public FigureOfMerit {
public:

   /**
    * \name Constructors with weights
    *
    * @{
    *
    * Constructor of a \f$d\f$-dimensional lattice rule with \f$n\f$ points and
    * weights \f$\gamma_j\f$ = `gamma[j]` for \f$j = 1, 2, …, d\f$. `prime`
    * indicates whether \f$n\f$ is a prime number (<tt>true</tt>) or not
    * (<tt>false</tt>).
    */
   Discrepancy (long n, bool prime, double gamma[], int d);
   Discrepancy (long n, bool prime, const LatCommon::ProductWeights & gamma, int d);
   Discrepancy (long n, bool prime, const LatCommon::OrderDependentWeights & gamma, int d);
   /**
    * @}
    */

   /**
    * Similar to the above constructor, except that the dimension of
    * `gamma` is `dimg`, which is smaller than \f$d\f$. This is useful for
    * projections when `gamma[j]` is 0 for \f$j = \mathtt{dimg} + 1, …,
    * d\f$.
    */
   Discrepancy (long n, bool prime, int d, double gamma[], int dimg);

   /**
    * \name Constructors with r > 1
    *
    * @{
    *
    * Constructor of a lattice rule with \f$n\f$ points as above, except that
    * now \f$r > 1\f$.
    */
   Discrepancy (long n, long r, bool prime, double gamma[], int d);
   Discrepancy (long n, long r, bool prime, const LatCommon::ProductWeights & gamma, int d);
   Discrepancy (long n, long r, bool prime, const LatCommon::OrderDependentWeights & gamma,
                   int d);
   /**
    * @}
    */

   /**
    * Destructor.
    */
   virtual ~Discrepancy();

   /**
    * Returns the number of points \f$n\f$.
    */
   virtual long getNumPoints();

   /**
    * Returns `true` if \f$n\f$ is prime, `false` otherwise.
    */
   bool isPrimeN();

   /**
    * Returns `true` if \f$n\f$ is a power of 2, `false` otherwise.
    */
   bool isPower2N();

   /**
    * Returns the dimension of the points.
    */
   int getDimension();

   /**
    * Print the weights.
    */
   const std::string toString() const;

   /**
    * Precomputes the Fourier factors
    * \f[
    *   C_k = \sideset{}’\sum_{-N/2<h\le N/2}\; \frac{e^{i2\pi hk/n}}{|h|}
    * \f]
    * for all \f$k = 0, 1,…, n-1\f$.
    */
   void setFourier1 (double *C, long n, long N);

   /**
    * Precomputes the Fourier factors
    * \f[
    *   C_k = \sideset{}’\sum_{h\in\mathbb Z}\; \frac{e^{i2\pi hk/n}}{|h|^{\alpha}}
    * \f]
    * for all \f$k = 0, 1,…, n-1\f$.
    */
   void setFourier (double *C, long n, int alpha);

   /**
    * Computes the factors \f$\sigma_k(m, \ell)\f$ for \f$1
    * \le\ell\le\mathtt{z.size() - 1}\f$, according to the formula
    * {@link REF__OrdepBound_eq_sigma1 (eq.sigma1)}.
    */
   void computeSigma (const std::vector<long> & z, double* C, long k) const;
protected:

   /**
    * \name Set the weights
    *
    * @{
    *
    * Sets the weights \f$\gamma_j\f$ = `gamma[j]` for \f$j = 0, 1, 2, …, d\f$.
    * Element `gamma[0]` is unused. Sometimes we may set
    * \f$\gamma_{\emptyset}=1\f$, corresponding to the empty set.
    */
   void setWeights (double gamma[], int d);
   void setWeights (const LatCommon::ProductWeights & gamma);
   void setWeights (const LatCommon::OrderDependentWeights & gamma);
   /**
    * @}
    */

   long m_n;                      // n
   bool m_prime;                  // true if n is prime
   bool m_power2;                 // true if n is a power of 2
   int m_dim;                     // Dimension of lattice

   int m_dimGam;                  // Dimension of m_Gamma
   long m_r;                      // r
   long m_N;                      // N = f(n, r, el)
   double *m_Gamma;               // Weights gamma[j]
   double *m_C;                   // Precomputed Fourier factors
   double *m_frac;                // The fractions m_frac[i] = i/n
   long *ytemp;                   // Work vector of dim = d + 1
   mutable double **m_Sigma;      // Factors sigma_k(m,j) for order-dependent bound
private:

/**
 * Performs some initialization operations.
 */
void init (long n, bool prime, int dim, long r, double gamma[], int dimg);
};

}
#endif
