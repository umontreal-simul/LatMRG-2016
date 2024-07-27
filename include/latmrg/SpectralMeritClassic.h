#ifndef SPECTRALMERIT_CLASSIC_H
#define SPECTRALMERIT_CLASSIC_H
#include "FigureOfMerit.h"
#include "latcommon/Weights.h"
#include "latcommon/Normalizer.h"
#include <string>


namespace LatMRG {

/**
 * This class implements the computation of the shortest vector in a lattice.
 *
 * <div class="LatSoft-bigskip"></div>
 */
class SpectralMeritClassic : public FigureOfMerit {
public:

/**
 * Constructor for a lattice or point set with \f$n\f$ points in dimension
 * `dim`.
 */
SpectralMeritClassic (long n, int dim);

   /**
    * Destructor.
    */
   ~SpectralMeritClassic();

   /**
    * Returns the (constructor) parameters of this object as a string.
    */
   virtual const std::string toString() const;

   /**
    * Computes the weighted normalized length of the shortest vector in
    * the dual lattice, and returns it with a negative sign, to allow for
    * minimization. The functions `setCoordLimits()` and `setNormalizer()`
    * must be called at least once before calling `compute` for the first
    * time.
    * \remark **Richard:** Pourquoi n’appelle-t-on pas ces fonctions dans
    * le constructeur?
    *
    *  Calling `setWeights()` is optional.
    */
   double compute (long a[], int dim);

   /**
    * Set the coordinate limits for each projection order. In the
    * documentation for the `TestProjections` class:
    *
    * - `coordLimits[0]` corresponds to \f$k+1\f$;<br>- `coordLimits[i]`
    * corresponds to \f$t_i\f$ for \f$i \geqslant1\f$.
    *
    * The elements must satisfy `coordLimits[i]`
    * \f$\geqslant\f$ `coordLimits[i+1]` for \f$i \geqslant1\f$. In more
    * details:
    *
    * - `coordLimits[0]` is the minimum projection order to be considered;
    * it is usually set to 2;<br>- `coordLimits[1]` is the maximum
    * projection order to be considered for successive dimensions;<br>-
    * `coordLimits[i]` for \f$i \geqslant2\f$ is the maximum coordinate
    * index to be considered for non-successive projections of order
    * \f$i\f$.
    */
   void setCoordLimits (int* coordLimits);

   /**
    * The normalized length of the shortest vector in the dual for a given
    * projection is divided by the weight of this projection. This makes
    * the projection look as if it had an even shorter shortest vector,
    * which makes it more likely that a lattice will be rejected because
    * of this projection.
    */
   void setWeights (const LatCommon::Weights* weights);

   /**
    * After the length of the shortest vector in the dual is computed, it
    * is normalized using this normalizer.
    */
   void setNormalizer (LatCommon::Normalizer* normalizer);

   /**
    * For debugging. Default value: `false`.
    */
   void setDebug (bool debug);

   /**
    * If stationary is `true`, the lattice is assumed to be
    * projection-stationary. This means that only projections involving
    * coordinate 1 will be considered. Default: `false`.
    */
   void setProjectionStationary (bool stationary);

   /**
    * The best values of the figure of merit computed for each projection
    * order is stored after each call to the `compute` function. This
    * allows for the algorithm that computes the weighted normalized
    * length of the shortest vector in the dual to abort as soon as it
    * reaches a value below these current best values. Calling
    * `resetLimits()` makes the instance forget about any previously
    * computed values. The `compute` function calls `resetLimits()` when
    * it is called for a different dimension than at the previous call.
    */
   void resetLimits();
protected:

   long m_n;                      // n
   int m_dim;                     // Dimension of lattice

   const LatCommon::Weights* m_weights;
   LatCommon::Normalizer* m_normalizer;
   int* m_coordLimits;

   bool m_debug;
   bool m_stationary;

   /**
    * The algorithm that computes the weighted normalized length of the shortest
    * vector in the dual aborts whenever it reaches a value below these lower
    * limits. Each value in the array corresponds to a single dimension.
    */
   double* m_meritLowerLimits;
   int m_last_dim;
};

}
#endif
