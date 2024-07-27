#ifndef SEARCHERCBC_H
#define SEARCHERCBC_H
#include "FigureOfMerit.h"
#include "Searcher.h"


namespace LatMRG {

/**
 * This class implements searches to find good rank-1 lattices with respect
 * to a given figure of merit, using component-by-component (CBC) searches,
 * random or exhaustive for each component. That is, one searches for the
 * best lattice by varying only one component at a time for each dimension.
 * Once the best component has been found for a given dimension, then this
 * value is fixed and we pass to the search for the next component.
 *
 * The figure of merit object in the `SearcherCBC` constructor *must* set the
 * number of points \f$n\f$, the maximal dimension \f$s\f$ of the lattice and
 * possibly, the weight factors \f$\gamma_j\f$. Then the search program will
 * examine different lattices with \f$n\f$, \f$s\f$ and \f$\gamma_j\f$ fixed
 * in order to find the best amongst those examined.
 *
 * <div class="LatSoft-bigskip"></div>
 */
class SearcherCBC : public Searcher {
public:

/**
 * The number of points \f$n\f$, the dimension \f$s\f$, and possibly the
 * \f$s\f$ weight factors \f$\gamma_j\f$ must be given in `merit`. The
 * \f$n\f$ points of each lattice will be generated in the search.
 */
SearcherCBC (const FigureOfMerit* merit);

   /**
    * Destructor.
    *
    * `double exhaust (int s, long n);`
    *
    * Exhaustive CBC search to find the lattice with the best (the
    * smallest) figure of merit in dimension \f$s\f$. The search runs over
    * all values \f$ a_j = 1, 2, …, (n-1)\f$ for a given dimension
    * \f$j\f$. Once the best lattice has been found in dimension \f$j\f$,
    * that coefficient \f$a_j\f$ is fixed and the search runs over all
    * values of \f$a_{j+1}\f$ for the next dimension. The first component
    * \f$a_1\f$ of the generator is always set to 1; \f$a_0\f$ is unused.
    * The method returns the best value of the figure of merit in
    * dimension \f$s\f$.
    *
    * `double exhaustPrime (int s, long n);`
    *
    * Similar to  {@link #exhaust(int) exhaust(s)}, except that only
    * values of \f$a_j\f$ *relatively prime* to \f$n\f$ are considered.
    *
    * `double random (int s, long n, int k);`
    *
    * Random CBC search to find the lattice with the best (the smallest)
    * figure of merit in dimension \f$s\f$. \f$k\f$ random values
    * \f$a_j\f$ are examined for each dimension \f$j\f$. Each random
    * component \f$a_j\f$ takes values over the integers \f$1, 2, …,
    * (n-1)\f$. The first component \f$a_1\f$ of the generator is always
    * set to 1; \f$a_0\f$ is unused. The method returns the best value of
    * the figure of merit found in dimension \f$s\f$.
    *
    * `double randomPrime (int s, long n, int k);`
    *
    * Similar to  {@link #random(int,int) random(s, k)}, except that only
    * values of \f$a_j\f$ *relatively prime* to \f$n\f$ are considered.
    */
   virtual ~SearcherCBC();

   /**
    * Returns the best value of the figure of merit found in the last
    * search, in each dimension up to \f$s\f$. The values returned are
    * \f$V_j\f$, for \f$j = 1, 2, …, s\f$. \f$V_0\f$ is unused.
    */
   const std::vector<double>& getBestVals();
protected:

   std::vector<double> m_bestVals; // best values of merit in dim = 1,2,...,s
   double exhaust (int s, long n, bool relPrime);
   double random (int s, long n, int k, bool relPrime);
};

}
#endif