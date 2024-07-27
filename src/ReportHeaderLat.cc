/*
ReportHeaderLat.h for ISO C++
version
authors: Hicham Wahbi
Frédérik Rozon
Richard Simard
*/

#include "latcommon/Const.h"
#include "latmrg/ReportHeaderLat.h"
#include <string>

using namespace LatCommon;


namespace LatMRG
{

ReportHeaderLat::ReportHeaderLat (Writer * writer, LatConfig * config,
                                  LatMRG::IntLattice * lattice): ReportHeader (writer)
{
   m_config = config;
   m_lattice = lattice;
}


void ReportHeaderLat::printHeader ()
{
   m_writer->newLine ();
   m_writer->writeString (
      "===========================================================================");
   m_writer->newLine ();
   m_writer->writeString ("Generator's components:");
   m_writer->newLine ();
   m_writer->writeString ("Number of components: ");
   m_writer->writeInt (m_config->J);
   m_writer->newLine ();

   // Écriture des informations sur les composantes

   m_writer->beginTabbedSection ();
   m_writer->addTab ();

   for (int i = 0; i < m_config->J; i++) {
      if (m_config->J > 1) {
         m_writer->newLine ();
         m_writer->writeString ("Component ");
         m_writer->writeInt (i + 1);
         m_writer->newLine ();

      }

      m_writer->writeString ("   Lattice type: ");
      m_writer->writeString (toStringGen (m_config->genType[i]));
      m_writer->newLine ();

      m_writer->writeString ("   m = ");
      m_writer->writeMScal (m_config->comp[i]->getM ());

      m_writer->newLine ();

      m_writer->writeString ("   k = ");
      m_writer->writeInt (m_config->comp[i]->k);
      m_writer->newLine ();

      if (m_config->comp[i]->k < 13) {
         // Write all coefficients of a
         for (int ii = 1; ii <= m_config->comp[i]->k; ii++) {
            m_writer->writeString ("      a_");
            m_writer->writeInt (ii);
            m_writer->writeString (" = ");
            m_writer->writeMScal (m_config->comp[i]->a[ii]);
            m_writer->newLine ();
         }
      } else {
         // Write only the non-zero coefficients of a
         for (int ii = 1; ii <= m_config->comp[i]->k; ii++) {
            if (0 != m_config->comp[i]->a[ii]) {
               m_writer->writeString ("      a_");
               m_writer->writeInt (ii);
               m_writer->writeString (" = ");
               m_writer->writeMScal (m_config->comp[i]->a[ii]);
               m_writer->newLine ();
            }
         }
         m_writer->writeString ("         All other a_j are  0");
         m_writer->newLine ();
      }

      if (m_config->J > 1) {
         m_writer->writeString ("   nj = ");
         m_writer->writeMScal (m_config->comp[i]->nj);
         m_writer->newLine ();
         m_writer->writeString ("   rho = ");
         m_writer->writeMScal (m_config->comp[i]->rho);
         m_writer->newLine ();
      }
   }

   if (m_config->J > 1) {
      m_writer->newLine ();
      m_writer->writeString ("Combined MRG:");
      m_writer->newLine ();
      m_writer->writeString ("   m = ");
      m_writer->writeMScal (m_lattice->getM ());
      m_writer->newLine ();
      m_writer->writeString ("   k = ");
      m_writer->writeInt (m_lattice->getOrder ());
      m_writer->newLine ();

      for (int i = 1; i <= m_lattice->getOrder (); i++) {
         m_writer->writeString ("      a_");
         m_writer->writeInt (i);
         m_writer->writeString (" = ");
         m_writer->writeMScal (m_lattice->getCoef ()[i]);
         m_writer->newLine ();
      }
      m_writer->writeString ("   rho = ");
      m_writer->writeMScal (m_lattice->getRho ());
      m_writer->newLine ();
      m_writer->writeString ("   lossRho = ");
      m_writer->writeMScal (m_lattice->getLossRho ());
      m_writer->newLine ();
   }

   m_writer->newLine ();
   m_writer->endTabbedSection ();

   m_writer->writeString ("Test: ");
   m_writer->writeString (toStringCriterion (m_config->criter));
   m_writer->newLine ();
   if (m_config->criter == PALPHA) {
      m_writer->writeString ("CalcType: ");
      m_writer->writeString (toStringCalc (m_config->calcPalpha));
      m_writer->newLine ();
   }

   if (m_config->criter == SPECTRAL) {
      m_writer->writeString ("Norm: ");
      m_writer->writeString (toStringNorm (m_config->norm));
      m_writer->newLine ();
      m_writer->writeString ("Normalisation: ");
      m_writer->writeString (toStringNorma (m_config->norma));
      m_writer->newLine ();
   }

   if (m_config->criter == PALPHA) {
      m_writer->writeString ("Prime m: ");
      m_writer->writeBool (m_config->primeM);
      m_writer->newLine ();
      m_writer->writeString ("Max period: ");
      m_writer->writeBool (m_config->maxPeriod);
      m_writer->newLine ();
      m_writer->writeString ("Alpha: ");
      m_writer->writeDouble (m_config->alpha);
      m_writer->newLine ();
      m_writer->writeString ("Seed: ");
      m_writer->writeInt (m_config->seed);
      m_writer->newLine ();
      m_writer->writeString ("Beta: { ");
      m_writer->writeDouble (m_config->Beta[0]);
      for (int i = 0; i <= m_config->td[1]; i++) {
         m_writer->writeString (", ");
         m_writer->writeDouble (m_config->Beta[i]);
      }
      m_writer->writeString (" }");

   } else {
      m_writer->writeString ("Lattice Type: ");
      m_writer->writeString (toStringLattice (m_config->latType));
      m_writer->newLine ();
      if (m_config->dualF)
         m_writer->writeString ("Lattice: DUAL");
      else
         m_writer->writeString ("Lattice: PRIMAL");
      m_writer->newLine ();
      m_writer->writeString ("Dimensions td:");
      for (int i = 0; i <= m_config->d; i++) {
         m_writer->writeString (" ");
         m_writer->writeInt (m_config->td[i]);
      }
      m_writer->newLine ();
   }

   if (m_config->lacunary) {
      int dim = m_config->td[1];
      m_writer->writeString ("dim = ");
      m_writer->writeInt (dim);
      m_writer->writeString ("\n\nLacunary = {   ");
      // Print indices by group of s
      MScal pre;
      conv (pre, -9);
      for (int i = 1; i <= dim; i++) {
         MScal r = m_config->Lac[i];
         if (pre < r - 1)
            m_writer->writeString ("\n   ");
         m_writer->writeMScal (r);
         m_writer->writeString ("    ");
         pre = r;
      }
      m_writer->writeString ("\n}\n");
   }

   m_writer->newLine ();
   m_writer->newLine ();
}

}
