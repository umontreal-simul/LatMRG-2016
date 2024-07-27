#include "latmrg/AverageCaseMerit.h"
#include "latmrg/ProjIteratorDefault.h"
#include <sstream>

//!#define DEBUG

#ifdef DEBUG
#include <iostream>
#endif

using namespace LatCommon;


namespace LatMRG
{

AverageCaseMerit::AverageCaseMerit (const ProjectionMerit & merit,
                         const Weights & weights, bool alwaysLastCoord)
      : m_merit (merit), m_weights (weights),
        m_alwaysLastCoord (alwaysLastCoord)
{}

//===========================================================================

const std::string AverageCaseMerit::toString () const
{
   std::ostringstream os;
   os << "average-case figure of merit based on " << m_merit.toString ()
   << " with " << m_weights;
   return os.str ();
}

//===========================================================================

double AverageCaseMerit::compute (const std::vector <long>&a, int n,
                                  double abortThreshold) const
{
   // FIXME: start at order 2
   int dim = a.size () - 1;
   ProjIteratorDefault proj (1, dim, 1, dim, false, m_alwaysLastCoord);
   return compute (a, n, proj, abortThreshold);
}

//===========================================================================

double AverageCaseMerit::compute (const std::vector <long>&a, int n,
                                  ProjIterator & proj, double abortThreshold) const
{
#ifdef DEBUG
   std::cout << "computing average-case merit for generator [";
   for (std::vector <long>::const_iterator it = a.begin () + 1;
         it != a.end (); ++it)
      std::cout << " " << *it;
   std::cout << " ]" << std::endl;
#endif

   int nproj = 0;
   double sum = 0.0;

   while (proj) {
      double weight = m_weights.getWeight (*proj);
      if (weight == 0.0) {
#ifdef DEBUG
         std::cout << "  skipping projection:";
         for (std::vector < int >::const_iterator it = proj->begin ();
               it != proj->end (); ++it)
            std::cout << " " << *it;
         std::cout << std::endl;
#endif
         ++proj;
         continue;
      }

#ifdef DEBUG
      std::cout << "  processing projection:";
      for (std::vector < int >::const_iterator it = proj->begin ();
            it != proj->end (); ++it)
         std::cout << " " << *it;
      std::cout << ":" << std::endl;
      std::cout << "    weight: " << weight << std::endl;
#endif

      double merit = m_merit.compute (a, n, *proj);
#ifdef DEBUG
      std::cout << "    merit:  " << merit << std::endl;
#endif

      sum += weight * merit;
      ++proj;

      if (abortThreshold > 0 && sum > abortThreshold) {
#ifdef DEBUG
         std::cout << "    threshold reached: aborting" << std::endl;
#endif
         break;
      }

      ++nproj;
   }

#ifdef DEBUG
   std::cout << "  projections considered: " << nproj << std::endl;
   std::cout << "  final merit: " << sum << std::endl;
#endif

   return sum;
}

//===========================================================================
}
