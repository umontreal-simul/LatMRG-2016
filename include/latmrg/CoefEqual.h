#ifndef COEFEQUAL_H
#define COEFEQUAL_H
#include "Coefficient.h"


namespace LatMRG {

/**
 * This class chooses the coefficients \f$a_i\f$ of a recursive generator of
 * order \f$k\f$ such that all the coefficients are equal by groups. There
 * are \f$s\f$ groups: the first group of \f$k_1\f$ coefficients are all
 * equal, the second group of \f$k_2 - k_1\f$ coefficients are all equal, and
 * so on until the last group of \f$(k - k_{s-1})\f$ coefficients which are
 * all equal. Thus the recursion is of the form
 * \f{align*}{
 *    x_n 
 *    & 
 *   =
 *  \alpha_1(x_{n-1} + \cdots+ x_{n-k_1}) + \alpha_2(x_{(n-1 - k_1)} + \cdots+ x_{n-k_2}) + 
 *  \\  &   
 *   \cdots+ \alpha_n (x_{(n-1 - k_{s-1})} + \cdots+ x_{n-k}) \bmod m.
 * \f}
 * For example, for a MRG of order \f$k=7\f$, one may consider only
 * recurrences of the form
 * \f[
 *   x_n = \alpha(x_{n-1} + x_{n-2}) + \beta(x_{n-3} + x_{n-4} + x_{n-5}) + \gamma(x_{n-6} + x_{n-7}) \bmod m,
 * \f]
 * where \f$\alpha, \beta, \gamma\f$ are arbitrary integers smaller than
 * \f$m\f$. An extreme case is when all coefficients are the same, as in
 * \f[
 *   x_n = \alpha(x_{n-1} + x_{n-2} + x_{n-3} + x_{n-4} + x_{n-5} + x_{n-6} + x_{n-7}) \bmod m.
 * \f]
 */
class CoefEqual: public Coefficient {
public:

/**
 * Constructor. The vector `I` (with \f$s+1\f$ elements), contains the
 * indices \f$k_0, k_1, k_2, …, k_s\f$ of the last member of each group of
 * equal coefficients, starting with `I[0]` \f$ = 0\f$ and ending with `I[s]`
 * \f$ = k\f$. For example, for a recurrence of order \f$k=7\f$ of the form
 * \f[
 *   x_n = \alpha(x_{n-1} + x_{n-2}) + \beta(x_{n-3} + x_{n-4} + x_{n-5}) + \gamma(x_{n-6} + x_{n-7}) \bmod m,
 * \f]
 * and thus with \f$s=3\f$ groups of equal coefficients, `I` must be \f$(0,
 * 2, 5, 7)\f$.
 */
CoefEqual (int *I, int s);

   /**
    * Destructor.
    */
   ~CoefEqual();

   /**
    * Sets the \f$r\f$-<em>th</em> group of equal coefficients `A[j]` \f$
    * = a_j\f$ to the value \f$q\f$, for \f$j = i, i-1, …, (k_{r-1}+
    * 1)\f$. The input value of \f$i=k_r\f$, the index of the upper member
    * of the \f$r\f$-<em>th</em> group. On return, the value of \f$i\f$ is
    * reset to \f$i = k_{r-1} + 1\f$, the index of the lowest member of
    * the \f$r\f$-<em>th</em> group.
    */
   void set (const MScal & q, MVect & A, int & i);
private:

/**
 * Contains the indices as defined in the constructor. It is a pointer to the
 * vector `I` defined elsewhere.
 */
int *m_I;

   /**
    * The number of different groups of unequal coefficients (it is equal
    * to the value of \f$s\f$ in the constructor).
    */
   int m_s;
};

}
#endif
