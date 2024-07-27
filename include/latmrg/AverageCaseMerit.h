#ifndef AVERAGE_CASE_MERIT_H
#define AVERAGE_CASE_MERIT_H
#include "latmrg/FigureOfMerit.h"
#include "latmrg/ProjectionMerit.h"
#include "latmrg/ProjIterator.h"
#include "latcommon/Weights.h"
#include <vector>


namespace LatMRG {

/**
 * This class implements generic average-case figures of merit of the form
 * \f[
 *   D^2(P_n) = \sum_{\mathfrak u\in\{1,…,s\}} \gamma_{\mathfrak u} D_{\mathfrak u}^2(P_n(\mathfrak u)),
 * \f]
 * for given weights \f$\gamma_{\mathfrak u}\f$ and using a per-projection
 * figure of merit \f$D_{\mathfrak u}^2\f$, for all projections
 * \f$\mathfrak u\f$.
 *
 */
class AverageCaseMerit : public FigureOfMerit {
public:

/**
 * Constructor. If `alwaysLastCoord` is `true`, a variant of the figure of
 * merit is created where only projections that include the last coordinate
 * are considered (this allows for faster CBC searches).
 */
AverageCaseMerit (const ProjectionMerit & merit, const LatCommon::Weights & weights,
                     bool alwaysLastCoord = false);

   /**
    * Destructor.
    */
   virtual ~AverageCaseMerit() {}

   /**
    * Returns a string representation of the figure of merit.
    */
   virtual const std::string toString() const;
virtual double compute (const std::vector<long>& a, int n,
                           double abortThreshold = 0) const;

/**
 * See `FigureOfMerit`. The implementation assumes that the internal
 * `ProjectionMerit` instance returns non-negative values only.
 */
virtual double compute (const std::vector<long>& a, int n,
                           ProjIterator & proj, double abortThreshold = 0) const;
protected:

   const ProjectionMerit & m_merit;
   const LatCommon::Weights & m_weights;
   const bool m_alwaysLastCoord;
};

}
#endif
