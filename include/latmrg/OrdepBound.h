#ifndef ORDEPBOUND_H
#define ORDEPBOUND_H
#include "BoundJSstar.h"


namespace LatMRG {

/**
 * This class computes bounds on the weighted star discrepancy for lattice
 * rules constructed with the CBC algorithm as
 * \anchor REF__OrdepBound_ordepbound proposed by Joe and Sinescu
 * \cite rSIN07a, \cite rSIN08b&thinsp;. Here we consider general weights of
 * a particular type, the so-called order-dependent weights which are
 * implemented in this class. For a lattice rule with \f$n\f$ points, the
 * bound on the weighted star discrepancy for general weights is given by
 * \f{align*}{
 *    e_{n,d}^2(\mathbf{z}) 
 *    & 
 *   =
 *  \frac{1}{n} \sum_{\mathbf{u} \subseteq\{1,2,…,d\}} \boldsymbol {\gamma}_{\mathbf{u}}\left(\sum_{k=0}^{n-1} \prod_{j\in\mathbf{u}} \left( \sideset{}’\sum_{-N/2<h\le N/2}\frac{e^{i2\pi hkz_j/n}}{|h|}\right)\right), 
 *  \\ 
 *    D^*_{n,\boldsymbol {\gamma}}(\mathbf{z}) 
 *    & 
 *   =
 *  \frac{e_{n,d}^2(\mathbf{z})}{2} + \frac{1}{N} \max_{\mathbf{u} \subseteq\{1,2,…,d\}}| \mathbf{u}|\boldsymbol {\gamma}_{\mathbf{u}}.
 * \f}
 * where \f$i=\sqrt{-1}\f$, \f$N = nr\f$, and \f$\mathbf{z}\f$ is the
 * generating vector of the lattice rule. Each component \f$z_j\f$ of the
 * generating vector is taken from the set \f$\{1, 2, …, n-1\}\f$ such that
 * \f$z_j\f$ is relatively prime with \f$n\f$. By
 * \f$\boldsymbol {\gamma}\f$ we denote the weights. The symbol
 * \f$\sum^{\prime}\f$ means that the \f$h=0\f$ term is omitted from the
 * sum.
 *
 * In practice, order-dependent weights are important since they allow a
 * significant reduction of the computational costs. Ordep-dependent weights
 * are weights that depend only on the number of elements of each subset of
 * \f$\{1,2,…,d\}\f$ so that subsets having the same cardinality have equal
 * weights. These weights are denoted by
 * \f$\Gamma_1,\Gamma_2,…,\Gamma_d\f$, where by \f$\Gamma_{\ell}\f$ we
 * mean the weight associated with any subset having \f$\ell\f$ elements.
 * For order-dependent weights, the bound is given by
 * \anchor REF__OrdepBound_eq_sinescu2
 * \f{align}{
 *    e_{n,d}^2(\mathbf{z}) 
 *    & 
 *   =
 *  \frac{1}{n} \sum_{\ell=1}^d \Gamma_{\ell}\sum_{\substack {\mathbf{u} \subseteq\{1,2,…,d\} \\
 *   |\mathbf{u}|=\ell
 *   }} \; \sum_{k=0}^{n-1}\prod_{j\in\mathbf{u}}\left( \sideset{}’\sum_{-N/2<h\le N/2}\; \frac{e^{i2\pi hkz_j/n}}{|h|}\right), \tag{eq.sinescu2} 
 *  \\ 
 *    D^*_{n,\boldsymbol {\gamma}}(\mathbf{z}) 
 *    & 
 *   =
 *  \frac{e_{n,d}^2(\mathbf{z})}{2} + \frac{1}{N} \max_{\ell\in\{1, 2,…,d\}}\ell\Gamma_{\ell}. \nonumber
 * \f}
 * In order to compute equation {@link REF__OrdepBound_eq_sinescu2
 * (eq.sinescu2)}, we see that we can write
 * \f[
 *   e_{n,d}^2(\mathbf{z})=\frac{1}{n}\sum_{k=0}^{n-1}\sum_{\ell=1}^d\Gamma_{\ell}\sigma_k(d, \ell),
 * \f]
 * where
 * \anchor REF__OrdepBound_eq_sigma1
 * \f[
 * \sigma_k(d,\ell)=\sum_{\substack {\mathbf{u} \subseteq\{1,2,…,d\} \\| \mathbf{u}|=\ell}}\prod_{j\in\mathbf{u}}C_k(z_j)\quad for\quad1\le\ell\le d, \tag{eq.sigma1}
 * \f]
 * with
 * \f[
 *   C_k(z)=\sideset{}’\sum_{-N/2<h\le N/2}\frac{e^{i2\pi hkz/n}}{|h|}.
 * \f]
 * Then we have a recursive formula to compute the quantities
 * \f$\sigma_k(d,\ell)\f$, that is
 * \f[
 * \sigma_k(d,\ell)=\sigma_k(d-1,\ell)+C_k(z_d)\sigma_k(d-1,\ell-1),
 * \f]
 * for \f$\ell=2,3,…,d-1\f$, while
 * \f$\sigma_k(d,1)=\sum\limits_{j=1}^dC_k(z_j)\f$ and
 * \f$\sigma_k(d,d)=\prod\limits_{j=1}^dC_k(z_j)\f$. 
 *
 * <div class="LatSoft-bigskip"></div>
 */
class OrdepBound: public BoundJSstar {
public:

/**
 * Constructor. The bound {@link REF__OrdepBound_eq_sinescu2 (eq.sinescu2)}
 * will be computed for rank-1 lattice rules with \f$n\f$ points and weights
 * \f$\Gamma_j\f$ = `gamma[j]` for \f$j = 1, 2, …, d\f$. `prime` indicates
 * whether \f$n\f$ is a prime number (<tt>true</tt>) or not (<tt>false</tt>).
 * In this case, \f$N=n\f$ in eq. {@link REF__OrdepBound_eq_sinescu2
 * (eq.sinescu2)}.
 */
OrdepBound (long n, bool prime, double gamma[], int d);

   /**
    * Similar to the above constructor, except that the dimension of
    * `gamma` is `dimg`, which is smaller than \f$d\f$. This is useful for
    * projections when `gamma[j]` is 0 for \f$j = \mathtt{dimg} + 1, …,
    * d\f$.
    */
   OrdepBound (long n, bool prime, int d, double gamma[], int dimg);

   /**
    * Constructor for a lattice rule with \f$n\f$ points as above, except
    * that now \f$nr=N\f$.
    */
   OrdepBound (long n, long r, bool prime, double gamma[], int d);

   /**
    * Destructor.
    */
   ~OrdepBound();

   /**
    * Returns the (constructor) parameters of this object as a string.
    */
   const std::string toString() const;

   /**
    * Computes and returns the bound {@link REF__OrdepBound_eq_sinescu2
    * (eq.sinescu2)} in dimension \f$d\f$ for the lattice with generating
    * vector \f$a_j =\f$ `a[j]`. Element `a[0]` is unused. Parameter `n`
    * must be the same value as used during object initialization
    * \remark **David:**
    *
    * this must be fixed. Parameter `abortThreshold` is ignored.
    */
   double compute (const std::vector<long> & a, int n, double abortThreshold) const;
private:

/**
 * Initialization. Computes the Fourier factors \f$C_k\f$.
 */
void init ();
};

}
#endif // ORDERBOUND_H