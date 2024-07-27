#ifndef POLY_H
#define POLY_H
#include "latcommon/Types.h" 
#include "IntFactorization.h"
#include "IntPrimitivity.h"
#include <string>


namespace LatMRG {

/**
 * This class implements polynomials \f$P(x)\f$ in \f$\mathbb Z_m[X]\f$
 * defined as
 * \anchor REF__PolyPE_eq_poly1
 * \f[
 *   P(x) = c_0 + c_1x^1 + c_2 x^2 + \cdots+ c_nx^n \tag{eq.poly1}
 * \f]
 * with degree \f$n\f$ and integer coefficients \f$c_i\f$ in
 * \f$\mathbb Z_m\f$. The arithmetic operations on objects of this class are
 * done modulo \f$m\f$ and modulo a polynomial \f$f(x)\f$ of degree \f$k\f$.
 * Thus all polynomials will be reduced modulo \f$f(x)\f$. In LatMRG, the
 * modulus polynomial \f$f(x)\f$ is usually written in the form
 * \anchor REF__PolyPE_eq_poly2
 * \f[
 *   f(x) = x^k - a_1x^{k-1} - \cdots- a_{k-1} x - a_k, \tag{eq.poly2}
 * \f]
 * and is associated with the recurrence
 * \anchor REF__PolyPE_eq_rec2
 * \f[
 *   x_n = (a_1x_{n-1} + a_2x_{n-2} + \cdots+ a_k x_{n-k}) \bmod m. \tag{eq.rec2}
 * \f]
 * The two functions `setM` and `setF` *must* be called to initialize the
 * modulus \f$m\f$ and the modulus polynomial \f$f(x)\f$ before doing any
 * arithmetic operations on `PolyPE` objects, otherwise the results are
 * unpredictable.
 *
 * Type `MScal` is used to represent polynomial coefficients. It may be
 * chosen as `long` for \f$m < 2^{50}\f$ (on 64-bit machines), or as the big
 * integer type `ZZ` otherwise. The possible associated types `MVect` are
 * `long*` and <tt>vec_ZZ</tt>. Type `PolE` for the polynomials may be chosen
 * as <tt>zz_pE</tt> when \f$m < 2^{50}\f$, or it may be set to
 * <tt>ZZ_pE</tt> which is implemented with the big integer type
 * <tt>ZZ_p</tt>.
 *
 */
class PolyPE : public PolE {
public:

/**
 * Initializes the modulus \f$m\f$ for this class. This must be called before
 * doing any operations on `PolyPE` objects, otherwise the results are
 * unpredictable.
 */
static void setM (const MScal & m);

   /**
    * Returns a read-only reference to \f$m\f$.
    */
   static const MScal & getM () { return m_m; }

   /**
    * Initializes the modulus polynomial \f$f(x) = c_0 + c_1x +c_2x^2 +
    * \cdots+ c_kx^k\f$ of degree \f$k\f$ for this class from the
    * coefficients \f$c_i = \f$<tt>C[i]</tt> of vector `C` of dimension
    * \f$k + 1\f$. This function must be called before doing any
    * arithmetic operations on `PolyPE` objects, otherwise the results are
    * unpredictable.
    */
   static void setF (const MVectP & C);

   /**
    * Same as above, but instead from a `MVect`.
    */
   static void setF (const MVect & C);

   /**
    * Returns the modulus polynomial \f$f(x)\f$.
    */
   static const PolX & getF () { return modulus().val(); }

   /**
    * Returns the degree of the modulus \f$f(x)\f$.
    */
   static long getK () { return degree(); }

   /**
    * Given a vector \f$C = [c_0, c_1, …, c_{k-1}, c_k]\f$, this function
    * reorders the components of \f$C\f$ in the form \f$C = [c_k, c_{k-1},
    * …, c_1, c_0]\f$ for `kind = 1`, and in the form \f$C =
    * [-c_k, -c_{k-1}, …, -c_1, 1]\f$ for `kind = 2`. For other values of
    * `kind`, it has no effect.
    */
   static void reverse (MVect & c, long k, int kind);

   /**
    * Minimal constructor: this object is set to the **0** polynomial.
    */
   PolyPE ();
   const PolX & getVal () { return rep(*this); }
   void setVal (long j);

   /**
    * Initializes this object to \f$C\f$.
    */
   void setVal (const MVect & C);

   /**
    * Initializes this object to the polynomial in `str`.
    */
   void setVal (std::string & str);

   /**
    * Sets \f$v = x^j \mod f(x) &nbsp;(\bmod&nbsp;m)\f$.
    */
   void powerMod (const MScal & j);

   /**
    * Returns the coefficients of this polynomial as a vector \f$C\f$ of
    * \f$k\f$ components, where \f$k\f$ is the degree of the modulus
    * \f$f(x)\f$.
    */
   void toVector (MVect & c);

   /**
    * Returns `true` if the modulus \f$f(x)\f$ is a primitive polynomial
    * modulo \f$m\f$. For this to be true, assuming that \f$f(x)\f$ has
    * the form {@link REF__PolyPE_eq_poly2 (eq.poly2)} above, the three
    * following conditions must be satisfied: \anchor REF__PolyPE_isprimi
    * <ol> <dt>None</dt>
    * <dd>
    * \f$[(-1)^{k+1} a_k]^{(m-1)/q} \bmod m \neq1\f$ for each prime
    * factor \f$q\f$ of \f$m - 1\f$;
    * </dd>
    * <dt>None</dt>
    * <dd>
    * \f$x^r \bmod(f(x),m) =  (-1)^{k+1} a_k \bmod m\f$; 
    * </dd>
    * <dt>None</dt>
    * <dd>
    * \f$x^{r/q} \bmod(f(x), m) \f$ has positive degree for each prime
    * factor \f$q\f$ of \f$r\f$, with \f$1<q< r\f$;
    * </dd>
    * </ol> where \f$r = (m^k - 1)/(m - 1)\f$. The factorizations of
    * \f$m-1\f$ and \f$r\f$ must be in `fm` and `fr` respectively.
    * Condition 1 is the same as saying that \f$(-1)^{k+1} a_k\f$ is a
    * primitive root of \f$m\f$. Condition 3 is automatically satisfied
    * when \f$r\f$ is prime.
    */
   bool isPrimitive (const IntPrimitivity & fm, const IntFactorization & fr);

   /**
    * Given the factorization of \f$r\f$, this method returns `true` if
    * conditions 2 and 3 above are satisfied by the modulus \f$f(x)\f$. It
    * does not check condition 1, assuming it to be `true`.
    */
   bool isPrimitive (const IntFactorization & r);

   /**
    * Returns this object as a string.
    */
   std::string toString () const;
private:

/**
 * Modulus of congruence.
 */
static MScal m_m;

   /**
    * Degree of the modulus polynomial \f$f\f$.
    */
   static long m_k;

   /**
    * The polynomial \f$f(x) = x\f$.
    */
   static PolX m_x;
};

}
#endif
