#ifndef PALPHA_H
#define PALPHA_H
#include "ProjectionMerit.h"


namespace LatMRG {

/**
 * This class implements the basic functions required to compute the
 * \f$P_{\alpha}\f$ figure of merit (see section  {@link
 * REF__sec1_sec_palpha palpha} ) with general weights for a lattice point
 * set \f$\Psi_s\f$ which is the intersection of a lattice \f$L\f$ and the
 * unit hypercube \f$[0, 1)^s\f$ in \f$s\f$ dimensions.
 * \f$\Psi_s\f$ contains \f$n\f$ points. For an arbitrary integer
 * \f$\alpha> 1\f$, it is defined as \cite mSLO94a, \cite vHIC98c,
 * \cite vHIC01a&thinsp;
 * \f[
 *   P_{2\alpha} (\Psi_s) = \sum_{\mathbf{0\neq h}\in L_s^*} \|\mathbf{h}\|^{-2\alpha},
 * \f]
 * where \f$L_s^*\f$ is the lattice dual to \f$L_s\f$, and the norm is
 * defined as \f$\|\mathbf{h}\| = \prod^s_{j=1} \max\{1, |h_j|\}\f$. When
 * \f$\alpha\f$ is integer, \f$P_{2\alpha}\f$ can be evaluated explicitly
 * as
 * \anchor REF__Palpha_palpha_1
 * \f[
 *   P_{2\alpha}(\Psi_s) = -1  +  \frac{1}{n}\sum_{\mathbf{u}\in\Psi_s}  \prod_{j=1}^s \left[1 - \frac{(-4\pi^2)^{\alpha}}{(2\alpha)!} B_{2\alpha}(u_j)\right], \tag{palpha.1}
 * \f]
 * where \f$B_{2\alpha}(x)\f$ is the Bernoulli polynomial of degree
 * \f$2\alpha\f$. The first Bernoulli polynomials of even degree are:
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
 * One may generalize the \f$P_{2\alpha}\f$ by introducing a weight for each
 * projection \cite vHIC98c&thinsp; (replace with reference with general
 * weights) to give the *weighted* \f$P_{2\alpha}\f$ defined by
 * \f[
 *   P_{2\alpha} (\Psi_s) = \sum_{\mathbf{0\neq h}\in L_s^*} \gamma_{\mathfrak u(\mathbf{h})} ^2 \|\mathbf{h}\|^{-2\alpha},
 * \f]
 * where \f$\gamma_{\mathfrak u}\f$ is the weight from projection
 * \f$\mathfrak u\f$ and
 * \f[
 * \mathfrak u(\vec h) = \left\{ j=1,…,s : h_j \neq0 \right\}.
 * \f]
 * For integer \f$\alpha\f$,
 * \f[
 *   P_{2\alpha}(\Psi_s) = \sum_{\emptyset\neq\mathfrak u\subseteq\{1,…,s\}} \gamma_{\mathfrak u} \; \frac{1}{n} \sum_{\mathbf{u} \in\Psi_s} \prod_{j \in\mathfrak u} \frac{-(-4\pi^2)^{\alpha}}{(2\alpha)!} B_{2\alpha}(u_j)
 * \f]
 * One recovers the original (unweighted) criterion {@link
 * REF__Palpha_palpha_1 (palpha.1)} for \f$P_{2\alpha}\f$ by setting
 * \f$\gamma_{\mathfrak u} = 1\f$ for all \f$\mathfrak u\in\{1,…,s\}\f$. 
 *
 * \remark This class uses functions from MyLib, part of TestU01.
 */
class Palpha : public ProjectionMerit {
public:

/**
 * Constructor. The \f$P_{\alpha}\f$ discrepancy {@link
 * REF__PalphaProduct_palpha_2 (palpha.2)}, with \f$\alpha=\f$ `alpha`, will
 * be computed for rank-1 lattices.
 */
Palpha (int alpha);

   /**
    * Destructor.
    */
   virtual ~Palpha() {}

   /**
    * Returns a string representation of the figure of merit.
    */
   virtual const std::string toString() const;

   /**
    * See class `ProjectionMerit`.
    */
   virtual double compute (const std::vector<long>& a, int n,
                           const std::vector<int>& projection) const;
protected:
   int m_alpha;                   // alpha

   class BernCache {
      // cache for scaled values of Bernoulli polynomials
      public:
    	    BernCache() { m_n = 0; m_alpha = 0; }
    	    void init(int alpha, int n); 
    	    double operator[](int i) const { return m_b[i % m_n]; }

      private:
    	    std::vector<double> m_b;
    	    int m_alpha;
    	    int m_n;
   };

   mutable BernCache m_bern;
};

}
#endif
