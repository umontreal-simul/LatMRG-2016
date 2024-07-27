#include "latmrg/WorstCaseMerit.h"
#include "latmrg/ProjIteratorDefault.h"
#include <sstream>
#include <cfloat>

// #define DEBUG

#ifdef DEBUG
#include <iostream>
#endif

using namespace LatCommon;


namespace LatMRG
{

WorstCaseMerit::WorstCaseMerit (const ProjectionMerit & merit,
                                const Weights & weights, bool alwaysLastCoord)
      : m_merit (merit), m_weights (weights),
        m_alwaysLastCoord (alwaysLastCoord)
{}

//===========================================================================

const std::string WorstCaseMerit::toString () const
{
   std::ostringstream os;
   os << "worst-case figure of merit based on " << m_merit.toString ()
   << " with " << m_weights;
   return os.str ();
}

//===========================================================================

double WorstCaseMerit::compute (const std::vector <long>&a, int n,
                                double abortThreshold) const
{
   // FIXME: start at order 2
   int dim = a.size () - 1;
   ProjIteratorDefault proj (1, dim, 1, dim, false, m_alwaysLastCoord);
   return compute (a, n, proj, abortThreshold);
}

//===========================================================================

double WorstCaseMerit::compute (const std::vector <long>&a, int n,
                                ProjIterator & proj, double abortThreshold) const
{
#ifdef DEBUG
   std::cout << "computing average-case merit for generator [";
   for (std::vector <long >::const_iterator it = a.begin () + 1;
         it != a.end (); ++it)
      std::cout << " " << *it;
   std::cout << " ]" << std::endl;
#endif

   double max = -DBL_MAX;

   while (proj) {

      double weight = m_weights.getWeight (*proj);
      if (weight == 0.0) {
#ifdef DEBUG
         std::cout << "  skipping projection:";
         for (std::vector <int>::const_iterator it = proj->begin ();
               it != proj->end (); ++it)
            std::cout << " " << *it;
         std::cout << std::endl;
#endif

         ++proj;
         continue;
      }
#ifdef DEBUG
      std::cout << "  processing projection:";
      for (std::vector <int>::const_iterator it = proj->begin ();
            it != proj->end (); ++it)
         std::cout << " " << *it;
      std::cout << ":" << std::endl;
      std::cout << "    weight: " << weight << std::endl;
#endif

      double merit = m_merit.compute (a, n, *proj);
#ifdef DEBUG

      std::cout << "    merit:  " << merit << std::endl;
#endif

      max = std::max (max, weight * merit);

      if (abortThreshold && max > abortThreshold) {
#ifdef DEBUG
         std::
         cout <<
         "  abort threshold reached: aborting computation of merit" <<
         std::endl;
#endif

         break;
      }

      ++proj;
   }

#ifdef DEBUG
   std::cout << "  final merit: " << max << std::endl;
#endif

   return max;
}

//===========================================================================
}
