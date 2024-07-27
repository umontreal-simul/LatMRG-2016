#ifndef PALPHA_PRODUCT_H
#define PALPHA_PRODUCT_H
#include "Discrepancy.h"
#include "latcommon/ProductWeights.h"


namespace LatMRG {

/**
 * This class implements the basic functions required to compute the
 * \f$P_{\alpha}\f$ figure of merit (see section  {@link
 * REF__sec1_sec_palpha palpha} ) for a lattice point set \f$\Psi_s\f$ which
 * is the intersection of a lattice \f$L\f$ and the unit hypercube \f$[0,
 * 1)^s\f$ in \f$s\f$ dimensions. \f$\Psi_s\f$ contains \f$n\f$ points. For
 * an arbitrary integer \f$\alpha> 1\f$, it is defined as \cite mSLO94a,
 * \cite vHIC98c, \cite vHIC01a&thinsp;
 * \f[
 *   P_{\alpha}(s) = \sum_{\mathbf{0\neq h}\in L_s^*} \|\mathbf{h}\|^{-\alpha},
 * \f]
 * where \f$L_s^*\f$ is the lattice dual to \f$L_s\f$, and the norm is
 * defined as \f$\|\mathbf{h}\| = \prod^s_{j=1} \max\{1, |h_j|\}\f$. When
 * \f$\alpha\f$ is even, \f$P_{\alpha}\f$ can be evaluated explicitly as
 * \anchor REF__PalphaProduct_palpha_1
 * \f[
 * \qquad P_{\alpha}(s) = -1  +  \frac{1}{n}\sum_{\mathbf{u}\in\Psi_s}  \prod_{j=1}^s \left(1 - \frac{(-1)^{\alpha/2}(2\pi)^{\alpha}}{\alpha!} B_{\alpha}(u_j)\right), \tag{palpha.1}
 * \f]
 * where \f$B_{\alpha}(x)\f$ is the Bernoulli polynomial of degree
 * \f$\alpha\f$. The first Bernoulli polynomials of even degree are:
 * \f{align*}{
 *    B_0(x) 
 *    & 
 *   =
 *    1
 *  \\ 
 *   B_2(x) 
 *    & 
 *   =
 *    x^2-x+1/6 
 *  \\ 
 *   B_4(x) 
 *    & 
 *   =
 *    x^4-2x^3+x^2-1/30
 *  \\ 
 *   B_6(x) 
 *    & 
 *   =
 *    x^6-3x^5+5x^4/2-x^2/2+1/42
 *  \\ 
 *   B_8(x) 
 *    & 
 *   =
 *    x^8-4x^7+14x^6/3 - 7x^4/3 +2x^2/3-1/30.
 * \f}
 * One may generalize the \f$P_{\alpha}\f$ by introducing a weight for each
 * dimension \cite vHIC98c&thinsp; to give the *weighted*
 * \f$P_{\alpha}\f$ defined by
 * \f[
 *   P_{\alpha}(s) = \sum_{\mathbf{0\neq h}\in L_s^*}\beta_I^2 \|\mathbf{h}\|^{-\alpha},
 * \f]
 * where the weights are such that
 * \f$\beta_I =\beta_0\prod_{j=1}^s\beta_j\f$, and for even \f$\alpha\f$
 * \anchor REF__PalphaProduct_palpha_2
 * \f[
 * \qquad P_{\alpha}(s) = \beta_0 \left\{-1  +  \frac{1}{n}\sum_{\mathbf{u}\in\Psi_s}   \prod_{j=1}^s \left(1 - \frac{(-1)^{\alpha/2}(2\pi)^{\alpha}\beta_j}{\alpha!} B_{\alpha}(u_j)\right) \right\}. \tag{palpha.2}
 * \f]
 * One recovers the original criterion {@link REF__PalphaProduct_palpha_1
 * (palpha.1)} for \f$P_{\alpha}\f$ by choosing all \f$\beta_j = 1\f$.
 *
 * \remark This class uses functions from MyLib, part of TestU01.
 */
class PalphaProduct : public Discrepancy {
public:

   PalphaProduct (int alpha, long n, bool prime, double beta[], int s);

/**
 * Constructors. The \f$P_{\alpha}\f$ bound {@link
 * REF__PalphaProduct_palpha_2 (palpha.2)}, with \f$\alpha=\f$ `alpha`, will
 * be computed for rank-1 lattices of \f$n\f$ points, in dimensions up to
 * \f$s\f$, with weights \f$\beta_j\f$ = `beta[j]` for \f$j = 1, 2, …, s\f$.
 * The flag `prime` must be set `true` if \f$n\f$ is a prime number.
 * \f$\beta_0\f$ is set to 1 for now.
 */
PalphaProduct (int alpha, long n, bool prime,
                  const LatCommon::ProductWeights & weights, int s);

   /**
    * Destructor.
    */
   ~PalphaProduct();

   /**
    * Returns the (constructor) parameters of this object as a string.
    */
   const std::string toString() const;

   /**
    * Returns the bound {@link REF__PalphaProduct_palpha_2 (palpha.2)} in
    * dimension of the vector `a` for the lattice with generator \f$z_j
    * =\f$ `a[j]`. Element `a[0]` is unused. Parameter `n` must be the
    * same value as used during object initialization.
    * \remark **David:**
    *
    * this must be fixed. Parameter `abortThreshold` is ignored.
    */
   virtual double compute (const std::vector<long> & a, int n,
                           double abortThreshold = 0) const;
protected:

/**
 * Initialization: computes the modified weights \f$G_s\f$ and the Bernoulli
 * factors \f$C_j\f$.
 */
void init (int alpha);

   /**
    * Precomputes the values of the Bernoulli polynomials \f$ C_j =
    * B_{\alpha}(j/n) \f$, for all \f$j = 0, 1,…, n-1\f$.
    */
   void setBernoulli (int alpha, long n);

   /**
    * Precomputes the modified weight factors \f$G_j =
    * (-1)^{\alpha/2} (2\pi)^{\alpha}\beta_j/\alpha!\f$ for \f$j = 1,
    * 2, …, s\f$.
    */
   void setG (int alpha, double beta[], int s);
   int m_alpha;                   // alpha
   double *m_G;                   // modified weights for j = 1, 2,..., s
};

}
#endif
