#ifndef INTPRIMITIVITY_H
#define INTPRIMITIVITY_H
#include <stdexcept>
#include "latcommon/Types.h"
#include "IntFactorization.h"


namespace LatMRG {

/**
 * This class deals with primitive roots and primitive elements modulo an
 * integer. Let \f$a\f$, \f$e\f$ and \f$p\f$ be integers, with \f$p\f$ a
 * prime number. Assume also that \f$a\f$ and \f$m=p^e\f$ are relatively
 * prime. The smallest positive integer \f$\lambda(m)\f$ for which
 * \f$a^{\lambda}= 1 \pmod\f$ is called the order of \f$a\f$ modulo
 * \f$m\f$. Any \f$a\f$ which has the maximum possible order for a given
 * \f$m\f$ is called a *primitive root* modulo \f$m\f$. For the following
 * important cases, the value of the order for given \f$m\f$ is
 * \cite rKNU98a&thinsp;:
 * \f{align*}{
 *  \lambda(2^e) 
 *    & 
 *   =
 *    2^{e-2}, \qquad e > 2 
 *  \\ 
 *  \lambda(p^e) 
 *    & 
 *   =
 *    p^{e-1}(p - 1), \qquad p > 2.
 * \f}
 */
class IntPrimitivity {
public:

   IntPrimitivity ();

   /**
    * Constructor fixing the modulus of congruence as \f$m = p^e\f$. The
    * argument \f$f\f$ must contain the prime factor decomposition of
    * \f$p-1\f$ and its inverse factors.
    */
   IntPrimitivity (const IntFactorization & f, const MScal & p, long e = 1);

   /**
    * Returns `true` if \f$a\f$ is a primitive element modulo \f$p^e\f$.
    * This method uses the prime factor decomposition of \f$p-1\f$ and its
    * inverse factors.
    */
   bool isPrimitiveElement (const MScal & a) const throw(std::range_error);

   /**
    * Returns `true` if \f$(-1)^{k+1}V[k]\f$ is a primitive element modulo
    * \f$p^e\f$. This method uses the prime factor decomposition of
    * \f$p-1\f$ and its inverse factors.
    */
   bool isPrimitiveElement (const MVect &V, int k) const throw(std::range_error);

   /**
    * Sets the value of \f$p\f$, \f$e\f$ and \f$m = p^e\f$.
    */
   void setpe (const MScal & p, long e);

   /**
    * Gets the value of \f$p\f$.
    */
   MScal getP () { return m_p; }

   /**
    * Gets the value of \f$e\f$.
    */
   long getE () { return m_e; }

   /**
    * Sets the value of \f$f\f$.
    */
   void setF (const IntFactorization & f) { m_f = f; }

   /**
    * Gets the value of \f$f\f$.
    */
   IntFactorization getF () { return m_f; }

   /**
    * Returns this object as a string.
    */
   std::string toString () const;
private:

/**
 * Prime number \f$p\f$.
 */
MScal m_p;

   /**
    * Exponent used to compute the modulus \f$m = p^e\f$.
    */
   long m_e;

   /**
    * The modulus \f$m = p^e\f$.
    */
   MScal m_m;

   /**
    * Factorization of \f$p-1\f$.
    */
   IntFactorization m_f;
};

}
#endif
