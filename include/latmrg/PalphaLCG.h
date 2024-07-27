#ifndef PALPHALCG_H
#define PALPHALCG_H 
#include "LatConfig.h"


namespace LatMRG {

/**
 * \remark **Richard:** Cette classe est basée sur les anciens programmes
 * Modula-2 pour un LCG. Il faudra peut-être songer à, soit l’éliminer, soit
 * la modifier, parce que ses méthodes semblent inutilement compliquées. Elle
 * est remplacée par la nouvelle classe `Palpha`. À vérifier.
 *
 *  This class implements the basic functions required to compute the
 * \f$P_{\alpha}\f$ figure of merit (see section  {@link
 * REF__sec1_sec_palpha palpha} ) for a *Korobov* lattice point set
 * \f$\Psi_s\f$, which is the intersection of a lattice \f$L\f$ and the unit
 * hypercube \f$[0, 1)^s\f$ in \f$s\f$ dimensions. See the doc in class
 * `Palpha` on page (FIXME: page#).
 * \remark **Richard:** Attention aux facteurs \f$\beta\f$: pas la même
 * définition que dans le nouveau `Palpha`.
 *
 * \remark This class uses functions from TestU01.
 */
class PalphaLCG {
public:

/**
 * Constructor.
 */
PalphaLCG (const LatConfig & config);

   /**
    * Destructor.
    */
   ~PalphaLCG ();

   /**
    * Sets the value of the multiplier to \f$a\f$.
    */
   void setMult (long a);

   /**
    * Computes and returns \f$P_2(s)\f$ in dimension \f$s\f$ for this
    * generator:
    * \f[
    * \qquad P_2(s) = \beta_0 \left(-1 + \frac{1}{m}\sum_{i=0}^{m-2} \prod_{j=0}^{s-1} \left(1 + 2\pi^2 \beta_{j+1}^2 \left(u^2_{i+j} - u_{i+j} + \frac{1}{6}\right)\right) + \frac{1}{m}\prod_{j=0}^{s-1} \left(1 + \frac{2\pi^2\beta_{j+1}^2}{6}\right) \right)
    * \f]
    * where \f$u_{i+j}\f$ is the number in \f$[0,1)\f$ returned by the
    * generator at step \f$(i+j)\f$-<em>th</em>, and \f$m\f$ is the
    * congruence modulus. Restrictions: the generator must have maximal
    * period, and one must have \f$\beta_j=1\f$ for all \f$j\f$; if these
    * conditions are not satisfied, one must use the slower
    * `calcPalpha2PerNonMax` method instead.
    */
   double calcPalpha2 (int s);

   /**
    * Computes and returns \f$P_2(s)\f$ in dimension \f$s\f$. The
    * generator may or may not have a maximal period. One makes sure to
    * get all points by considering all multiples \f$i\f$ of vector
    * \f$((1,a, a^2,…,a^{s-1})\mod m) /m\f$, where \f$i = 0, 1, 2 , …,
    * (m-1)\f$. This vector corresponds to \f$(u_i,
    * u_{i+1},...,u_{i+s-1})\f$ when the generator has maximal period.
    * This method gives the same result as `calcPalpha2` when all
    * \f$\beta_j=1\f$, but ensures that one does not get stuck in a cycle
    * when the period is not maximal.
    */
   double calcPalpha2PerNonMax (int s);

   /**
    * Computes and returns \f$P_2(s)\f$ in dimension \f$s\f$ for this
    * generator. This is a slow method for verification.
    */
   double calcPalpha2Verif (int s);

   /**
    * Computes \f$P_2(s)\f$ in all dimensions from `minDim` to `maxDim`
    * for this generator, assuming that all weights \f$\beta_j\f$ are 1,
    * and that the generator has maximal period. Returns the values
    * \f$P_2(s)\f$ in array `P2[s]`.
    */
   void calcPalpha2 (int minDim, int maxDim, double P2[]);

   /**
    * Computes \f$P_2(s)\f$ in all dimensions from `minDim` to `maxDim`
    * for this generator, assuming that all weights \f$\beta_j = 1\f$.
    * Returns the values \f$P_2(s)\f$ in array `P2[s]`.
    */
   void calcPalpha2PerNonMax (int minDim, int maxDim, double P2[]);

   /**
    * Similar to `calcPalpha2` but for \f$P_4(s)\f$:
    * \f[
    * \qquad P_4(s) = \beta_0 \left(-1 + \frac{1}{m}\sum_{i=0}^{m-2} \prod_{j=0}^{s-1} \left(1 - \frac{2\pi^4 \beta_{j+1}^4}{3} \left(u_{i+j}^4 - 2u_{i+j}^3 + u_{i+j}^2 - \frac{1}{30}\right)\right) + \frac{1}{m}\prod_{j=0}^{s-1} \left(1 + \frac{2\pi^4\beta_{j+1}^4}{90}\right) \right)
    * \f]
    */
   double calcPalpha4 (int s);

   /**
    * Similar to `calcPalpha2PerNonMax` but for \f$P_4(s)\f$.
    */
   double calcPalpha4PerNonMax (int s);

   /**
    * Similar to `calcPalpha2` but for \f$P_6(s)\f$:
    * \f{align*}{
    *    P_6(s) 
    *    & 
    *   = 
    *  \beta_0 \left(-1 + \frac{1}{m}\sum_{i=0}^{m-2} \prod_{j=0}^{s-1} \left(1 + \frac{4\pi^6 \beta_{j+1}^6}{45} \left(u_{i+j}^6 - 3u_{i+j}^5 + \frac{5}{2} u_{i+j}^4 - \frac{1}{2}u_{i+j}^2 + \frac{1}{42}\right)\right)\right. + 
    *  \\  &   
    *  \qquad\qquad\left.\frac{1}{m}\prod_{j=0}^{s-1} \left(1 + \frac{2\pi^6\beta_{j+1}^6}{945}\right) \right)
    * \f}
    */
   double calcPalpha6 (int s);

   /**
    * Similar to `calcPalpha2PerNonMax` but for \f$P_6(s)\f$.
    */
   double calcPalpha6PerNonMax (int s);

   /**
    * Similar to `calcPalpha2` but for \f$P_8(s)\f$:
    * \f{align*}{
    *  \qquad P_8(s) 
    *    & 
    *   = 
    *  \beta_0 \left(-1 + \frac{1}{m}\sum_{i=0}^{m-2} \prod_{j=0}^{s-1} \left(1 - \frac{2\pi^8 \beta_{j+1}^8}{315} \left(u_{i+j}^8 - 4u_{i+j}^7 + \frac{14}{3}u_{i+j}^6 - \frac{7}{3} u_{i+j}^4 + \frac{2}{3}u_{i+j}^2 - \frac{1}{30}\right)\right)\right. 
    *  \\  &   
    *  \qquad\quad{} + \left.\frac{1}{m}\prod_{j=0}^{s-1} \left(1 + \frac{2\pi^8\beta_{j+1}^8}{9450}\right) \right)
    * \f}
    */
   double calcPalpha8 (int s);

   /**
    * Similar to `calcPalpha2PerNonMax` but for \f$P_8(s)\f$.
    */
   double calcPalpha8PerNonMax (int s);
private:

   LatCommon::CalcType m_calcType;
   LatCommon::GenType m_genType;             // Generator type
   int m_m;                       // Modulus of congruence
   bool m_prime;                  // true if m_m is prime
   int m_a;                       // Multiplier
   int m_minDim;                  // Minimal dimension
   int m_maxDim;                  // Maximal dimension
   int m_alpha;  
   int m_seed;  
   double *m_Beta;                // Beta factors
   double *TabU;                  // Work variable
   double *Fac;

/**
 * Properties of the generator and of the computation.
 */
// Work variable
};

}
#endif
