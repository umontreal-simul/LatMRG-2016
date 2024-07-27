#ifndef BOUNDJSINTER_H
#define BOUNDJSINTER_H
#include "Discrepancy.h"


namespace LatMRG {

/**
 * This class computes bounds on the weighted star discrepancy for
 * intermediate-rank lattice rules constructed with the CBC algorithm as
 * proposed by Joe and Sinescu \cite rSIN08b, \cite rSIN07b&thinsp;. Here we
 * take intermediate-rank lattice rules with \f$N=\ell^rn\f$ points, where
 * \f$n\f$ is prime and \f$\ell\f$ is fixed such that
 * \f$\gcd(\ell,n)=1\f$. We also take a fixed \f$r\f$ from the set
 * \f$\{0,1,…,d\}\f$. For the particular values \f$r=0\f$ or
 * \f$\ell=1\f$, these intermediate-rank lattice rules become the usual
 * rank-\f$1\f$ lattice rules. In the case of product weights, the bound is
 * given by
 * \anchor REF__BoundJSInter_eq_sinescu3
 * \f{align}{
 *    e_{n,d}^2(\mathbf{z}) 
 *    & 
 *   =
 *  \frac{1}{n} \sum_{k=0}^{n-1} \prod_{j=1}^d \Bigg\{ \beta_j + \tilde{\gamma}_j \sideset{}’\sum_{-{\tilde{N}}_j/2<h\le{\tilde{N}}_j/2}\; \frac{e^{i2\pi hk\tilde{z}_j/n}}{|h|}\Bigg\} - \prod_{j=1}^d \beta_j \tag{eq.sinescu3} 
 *  \\ 
 *    D^*_{n,\boldsymbol {\gamma}}(\mathbf{z}) 
 *    & 
 *   =
 *  \frac{e_{n,d}^2(\mathbf{z})}{2} + \prod_{j=1}^d\beta_j-\prod_{j=1}^d\left(\beta_j-\frac{\gamma_j}{N}\right), \nonumber
 * \f}
 * where \f$\beta_j = 1 + \gamma_j\f$ for all \f$j=1,…,d\f$,
 * \f$i=\sqrt{-1}\f$, and \f$\mathbf{z}\f$ is the generating vector of the
 * lattice rule. Each component \f$z_j\f$ of the generating vector is taken
 * from the set \f$\{1, 2, …, n-1\}\f$ such that \f$z_j\f$ is relatively
 * prime with \f$n\f$. By \f$\boldsymbol {\gamma}\f$ we denote the weights.
 * The symbol \f$\sum^{\prime}\f$ means that the \f$h=0\f$ term is omitted
 * from the sum. In the above, we also used the following notations:
 * \f[
 * \tilde{\gamma}_j = \left\{ \begin{array}{ll}
 *  \gamma_j/\ell, 
 *    & 
 *  \quad\text{ for } 1 \le j \le r, 
 *  \\ 
 *  \gamma_j, 
 *    & 
 *  \quad\text{ for } r+1 \le j \le d, 
 * \end{array} \right.
 * \f]
 * \f[
 * \tilde{N}_j = \left\{ \begin{array}{ll}
 *    N/\ell, 
 *    & 
 *  \quad\text{ for } 1 \le j \le r, 
 *  \\ 
 *    N, 
 *    & 
 *  \quad\text{ for } r+1 \le j \le d, 
 * \end{array} \right.
 * \f]
 * \f[
 * \tilde{z}_j = \left\{ \begin{array}{ll}
 *  \ell z_j, 
 *    & 
 *  \quad \text{ for } 1 \le j \le r, 
 *  \\ 
 *    z_j, 
 *    & 
 *  \quad\text{ for } r+1 \le j \le d, 
 * \end{array} \right.
 * \f]
 * <div class="LatSoft-bigskip"></div>
 */
class BoundJSInter:  public Discrepancy {
public:

/**
 * Constructor. The bound {@link REF__BoundJSInter_eq_sinescu3 (eq.sinescu3)}
 * will be computed for intermediate-rank lattice rules with
 * \f$N=\ell^rn\f$ points and weights \f$\gamma_j\f$ = `gamma[j]` for \f$j
 * = 1, 2, …, d\f$. Element `gamma[0]` is unused. We also need to take
 * \f$n\f$ prime and \f$\ell\f$ chosen such that \f$\gcd(\ell,n)=1\f$.
 */
BoundJSInter (long n, long r, long el, double gamma[], int d);

   /**
    * Destructor.
    */
   ~BoundJSInter();

   /**
    * Computes and returns the bound {@link REF__BoundJSInter_eq_sinescu3
    * (eq.sinescu3)} in dimension \f$d\f$ for the lattice rule with
    * generating vector \f$z_j =\f$ `z[j]`, for \f$j = 1, 2, …, d\f$.
    * Element `z[0]` is unused.
    */
   double compute (long z[], int d);

   /**
    * Returns the (constructor) parameters of this object as a string.
    */
   const std::string toString() const;

   /**
    * Returns the number of points \f$N=\ell^rn\f$.
    */
   long getNumPoints();
private:

/**
 * Sets the values of \f$N=\ell^rn\f$ and \f$\tilde{N}=N/\ell\f$.
 */
void setN();

   /**
    * Sets the factors \f$\beta_j\f$ = `1 + gamma[j]`, for \f$j = 1, 2,
    * …, d\f$. Element \f$\beta_0\f$ is unused.
    */
   void setBeta (double gamma[], int d);

   /**
    * Initialization. Sets the factors \f$\tilde{\gamma}_j\f$ and \f$
    * \tilde{N}_j \f$, where \f$\gamma_j = \f$ `gamma[j]` and
    * \f$N=\ell^r n\f$.
    */
   void init (long el, double gamma[], int d);
   long m_el;                     // el
   long m_NTil;                   // tilde{N}
   long *m_zTil;                  // tilde{z_j}
   double *m_CTil;                // tilde{C_k}
   double *m_GammaTil;            // tilde{gamma_j}
   double *m_Beta;                // beta_j
};

}
#endif
