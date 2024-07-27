#ifndef EXACTDISCSTAR_H
#define EXACTDISCSTAR_H
#include "Discrepancy.h"


namespace LatMRG {

/**
 * This class computes the exact star discrepancy for lattice rules with
 * \f$n\f$ points in dimension \f$d\f$. The star discrepancy of a point set
 * \f$P_n\subseteq[0,1]^d\f$ is defined by (see for instance
 * \cite rNIE92b&thinsp;)
 * \anchor REF__ExactDiscStar_eq_stardis
 * \f[
 * \tag{eq.stardis} D^*(P_n) = \sup_{{\mathbf{x}}\in[0,1]^d}\left|\frac{\left|[\boldsymbol{0},\mathbf{x})\cap P_n\right|}{n}-\prod_{j=1}^dx_j\right|.
 * \f]
 * There is no efficient algorithm to compute the star discrepancy. More
 * precisely, any algorithm would require at least \f$O(n^d)\f$ operations,
 * which is unpractical when \f$n\f$ and \f$d\f$ are large. For small
 * dimensions however, one may use the algorithm proposed by Bundschuh and
 * Zhu \cite rBUN93&thinsp;.
 *
 * For \f$d=2\f$, the star discrepancy of the point set \f$P_n
 * =\{(x_k,y_k), 1\le k\le n\}\f$ is given by
 * \anchor REF__ExactDiscStar_eq_stardisd2
 * \f[
 * \tag{eq.stardisd2} D^*(P_n) = \max_{0\le\ell\le n}  \max_{0\le k\le\ell}  \max\left(\left|\frac{k}{n} -x_{\ell}\xi_{\ell,k}\right|, \left|\frac{k}{n}-x_{\ell+1}\xi_{\ell,k+1}\right|\right),
 * \f]
 * where we put \f$(x_0, y_0) = (0, 0)\f$ and \f$(x_{n+1}, y_{n+1}) =(1,
 * 1)\f$, while the coordinates \f$y_i\f$, for \f$i=0, 1, …, \ell, n+1\f$,
 * are rearranged in increasing order and rewritten as
 * \f$0=\xi_{\ell,0}\le\xi_{\ell,1}\le\cdots\le\xi_{\ell,\ell}<\xi_{\ell,\ell+1} =1\f$.
 *
 * A similar but slightly more complicated algorithm can be written for
 * \f$d=3\f$. Let’s consider the point set
 * \f$P_n=\{(x_k,y_k,z_k), 1\le k\le n\}\f$ and put \f$(x_0, y_0, z_0) = (0,
 * 0, 0)\f$ and \f$(x_{n+1}, y_{n+1}) =(1, 1, 1)\f$. Then, the coordinates
 * \f$y_i\f$, for \f$i=0, 1, …, \ell, n+1\f$, are rearranged in increasing
 * order and rewritten as above, while the corresponding coordinates
 * \f$z_i\f$ for \f$i=0, 1, …, \ell, n+1\f$ are rewritten as
 * \f$0=\eta_{\ell,t,0}\le\eta_{\ell,t,1}\le\cdots\le\eta_{\ell,t,t}< \eta_{\ell,t,t+1} =1\f$
 * for every fixed \f$0\le\ell\le n\f$ and for any \f$t=0,1,…,\ell\f$. The
 * star discrepancy is then given by
 * \anchor REF__ExactDiscStar_eq_stardisd3
 * \f[
 * \tag{eq.stardisd3} D^*(P_n) = \max_{0\le\ell\le n}  \max_{0\le t\le\ell}  \max_{0\le k\le t}   \max\left(\left|\frac{k}{n}-x_{\ell}\xi_{\ell,t}\eta_{\ell,t,k}\right|, \left|\frac{k}{n}-x_{\ell+1}\xi_{\ell,t+1}\eta_{\ell,t,k+1}\right|\right).
 * \f]
 * The algorithm can be extended for higher dimensions in a similar fashion.
 *
 * <div class="LatSoft-bigskip"></div>
 */
class ExactDiscStar: public Discrepancy {
public:

/**
 * Constructor. The star-discrepancy {@link REF__ExactDiscStar_eq_stardis
 * (eq.stardis)} will be computed for rank-1 lattice rules with \f$n\f$
 * points for \f$j = 1, 2, …, d\f$. `prime` indicates whether \f$n\f$ is a
 * prime number (<tt>true</tt>) or not (<tt>false</tt>). *WARNING:* Only
 * dimensions \f$d=2\f$ and \f$d=3\f$ are implemented.
 */
ExactDiscStar (long n, bool prime, int d);

   /**
    * Destructor.
    */
   ~ExactDiscStar();

   /**
    * Returns the (constructor) parameters of this object as a string.
    */
   const std::string toString() const;

   /**
    * Computes and returns the exact star discrepancy in dimension \f$d\f$
    * for the lattice rule with generating vector \f$z_j =\f$ `z[j]`, for
    * \f$j = 1, 2, …, d\f$. Element `z[0]` is unused. *WARNING:* Only
    * dimensions \f$d=2\f$ and \f$d=3\f$ are implemented.
    */
   double compute (long z[], int d);
private:

   void init (int d);

   void insert2Dim (double y, int n);
   void insert_yz (double y, double z, int n);
   void insert_z (double z, int n);

   double compute2Dim (long z[]);   // compute when dim = 2
   double compute3Dim (long z[]);   // compute when dim = 3

   double *m_x;          // x-coordinates
   double *m_y;          // y-coordinates
   double *m_z;          // z-coordinates
   double *m_Ksi;
   double *m_Eta;
   double *m_temp;
};

}
#endif