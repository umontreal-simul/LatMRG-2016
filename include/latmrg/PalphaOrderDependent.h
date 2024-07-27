#ifndef PALPHA_ORDER_DEPENDENT_H
#define PALPHA_ORDER_DEPENDENT_H
#include "Discrepancy.h"
#include "latcommon/OrderDependentWeights.h"


namespace LatMRG {

/**
 * This class computes a weighted discrepancy criterion for lattice rules
 * constructed with the CBC algorithm. The weighted discrepancy considered
 * here is based on the well known figure of merit \f$P_{\alpha}\f$, while
 * the weights considered here are order-dependent.
 *
 * For a lattice rule with \f$n\f$ points, this weighted
 * \f$P_{\alpha}\f$ discrepancy for order-dependent weights is given by
 * \anchor REF__PalphaOrderDependent_eq_ordeppalpha
 * \f{align}{
 *    D_{n,d,\boldsymbol {\gamma}}(\mathbf{z}, \alpha) 
 *    & 
 *   =
 *  \frac{1}{n} \sum_{\ell=1}^d \Gamma_{\ell}\sum_{\substack {\mathbf{u} \subseteq\{1,2,…,d\} \\
 *   |\mathbf{u}|=\ell
 *   }} \; \sum_{k=0}^{n-1}\prod_{j\in\mathbf{u}}\left( \sideset{}’\sum_{h\in\mathbb Z}\; \frac{e^{i2\pi hkz_j/n}}{|h|^{\alpha}}\right)\tag{eq.ordeppalpha}.
 * \f}
 * The order-dependent weights are denoted by
 * \f$\Gamma_1,\Gamma_2,…,\Gamma_d\f$, where by \f$\Gamma_{\ell}\f$ we
 * mean the weight associated with any subset having \f$\ell\f$ elements
 * [see also the bound on the weighted star discrepancy for order-dependent
 * weights in class `OrdepBound` (on page (FIXME: page#) of this guide)].
 *
 * In order to compute equation {@link
 * REF__PalphaOrderDependent_eq_ordeppalpha (eq.ordeppalpha)}, we can use the
 * same algorithm as for the star discrepancy with order-dependent weights,
 * the difference being the quantities \f$C_k(z)\f$ which are here given by
 * \anchor REF__PalphaOrderDependent_eq_ckz
 * \f[
 * \tag{eq.ckz} C_k(z, \alpha)=\sideset{}’\sum_{h\in\mathbb Z}\frac{e^{i2\pi hkz/n}}{|h|^{\alpha}}.
 * \f]
 * If \f$\alpha\f$ is an even integer, then we may use the Bernoulli
 * polynomials to express the quantities \f$C_k(z)\f$. For any other
 * \f$\alpha>1\f$, we could only obtain approximate values by truncating the
 * sum in {@link REF__PalphaOrderDependent_eq_ckz (eq.ckz)} (<em>...yet to be
 * done...</em>).
 *
 * Apart from this difference, the algorithm will follow the same idea as for
 * the weighted star discrepancy with order-dependent weights. Below we
 * present a few details of the implementation.
 *
 * <div class="LatSoft-bigskip"></div>
 */
class PalphaOrderDependent: public Discrepancy {
public:

   PalphaOrderDependent (int alpha, long n, bool prime, double gamma[], int d);

   /**
    * Constructors.
    */
   PalphaOrderDependent (int alpha, long n, bool prime,
                            const LatCommon::OrderDependentWeights & weights, int d);

   /**
    * Similar to the above constructor, except that the dimension of
    * `gamma` is `dimg`, which is smaller than \f$d\f$. This is useful for
    * projections when `gamma[j]` is 0 for \f$j = \mathtt{dimg} + 1, …,
    * d\f$.
    */
   PalphaOrderDependent (int alpha, long n, bool prime, int d,
                         double gamma[], int dimg);

   /**
    * Destructor.
    */
   ~PalphaOrderDependent();

   /**
    * Returns the (constructor) parameters of this object as a string.
    */
   const std::string toString() const;

   /**
    * Computes the weighted discrepancy {@link
    * REF__PalphaOrderDependent_eq_ordeppalpha (eq.ordeppalpha)} for
    * order-dependent weights. Element `a[0]` is unused. Parameter `n`
    * must be the same value as used during object initialization
    * \remark **David:**
    *
    * this must be fixed. Parameter `abortThreshold` is ignored.
    */
   virtual double compute (const std::vector<long>& a, int n,
                           double abortThreshold = 0) const;
protected:

/**
 * Does the initialization.
 */
void init (int alpha);
int m_alpha;

/**
 * This is the value of \f$\alpha\f$.
 */
// alpha
};

}
#endif
