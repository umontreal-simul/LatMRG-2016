#ifndef COEFZERO_H
#define COEFZERO_H
#include "Coefficient.h"


namespace LatMRG {

/**
 * This class chooses the coefficients \f$a_i\f$ of a recursive generator of
 * order \f$k\f$ such that all the coefficients are zero except for a select
 * few. There are \f$s\f$ non-zero coefficients: they have indices \f$k_1,
 * k_2, …, k_s\f$. Thus the recursion is of the form
 * \f[
 *   x_n = \alpha_1 x_{n-k_1} + \alpha_2x_{n-k_2} + \cdots+ \alpha_s x_{n-k} \bmod m.
 * \f]
 * For example, for a MRG of order \f$k=101\f$, one may consider only
 * recurrences with three non-zero coefficients
 * \f[
 *   x_n = (\alpha x_{n-3} + \beta x_{n-51} + \gamma x_{n-101}) \bmod m,
 * \f]
 * where \f$\alpha, \beta, \gamma\f$ are arbitrary integers smaller than
 * \f$m\f$. This allows for a simpler description of the input files and a
 * faster search for such generators than examining separately all zero
 * coefficients.
 *
 */
class CoefZero: public Coefficient {
public:

/**
 * Constructor. The vector `I` (with \f$s+1\f$ elements), contains the
 * indices \f$k_1, k_2, …, k_s\f$ of the non-zero coefficients, starting with
 * `I[1]` \f$ = k_1\f$ and ending with `I[s]` \f$ = k\f$. For example, for a
 * recurrence of order \f$k=7\f$ of the form
 * \f[
 *   x_n = (\alpha x_{n-3} + \beta x_{n-51} + \gamma x_{n-101}) \bmod m,
 * \f]
 * and thus with \f$s=3\f$ non-zero coefficients, `I` must be \f$(0, 3, 51,
 * 101)\f$.
 */
CoefZero (int *I, int s);

   /**
    * Destructor.
    */
   ~CoefZero();

   /**
    * Sets the \f$r\f$-<em>th</em> non-zero coefficient `A[i]` \f$ =
    * a_i\f$ to the value \f$q\f$, for \f$i = k_r\f$. On return, the value
    * of \f$i\f$ will be reset to \f$i = k_{r-1} + 1\f$, where \f$i =
    * k_{r-1}\f$ is the index of the next non-zero coefficient.
    */
   void set (const MScal & q, MVect & A, int & i);
private:

/**
 * Contains the indices as defined in the constructor. It is a pointer to the
 * vector `I` defined elsewhere.
 */
int *m_I;

   /**
    * The number of non-zero coefficients (it is equal to the value of
    * \f$s\f$ in the constructor).
    */
   int m_s;
};

}
#endif
