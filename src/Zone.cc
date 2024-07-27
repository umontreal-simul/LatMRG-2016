#include "latmrg/Zone.h"
#include "latcommon/Util.h"
#include "latmrg/SeekConfig.h"

#include <iostream>
#include <cassert>

using namespace std;
using namespace LatCommon;


namespace LatMRG
{

MMat Zone::Frontiere;

// Si DivQ[i] = true, on cherche les a=m/q avec q entre les bornes.
const bool Zone::DivQ[] = {true, false, true}; // For ZONE1, ZONE2, ZONE3



Zone::Zone ()
{
   inf = 0;
   sup = 0;
   No = ZONE3;
   frac = 0.0;
   smallF = false;
   supMsH = 0;
   T1 = 0;
   T2 = 0;
   nextZone = 0;
}

//===========================================================================

Zone::~Zone ()
{
   if (nextZone != 0) {
      delete nextZone;
      nextZone = 0;
   }
/*
   Zone *Z = nextZone;
   while (Z != 0) {
      Zone *p;
      p = Z->nextZone;
      delete Z;
      Z = p;
   }
*/
}


//===========================================================================

void Zone::initFrontieres (const SeekConfig & config)
{
   CreateMatr (Frontiere, config.J, NZONES);
   // Calcul des frontieres superieures des 3 zones
   for (int s = 0; s < config.J; s++) {
      Frontiere[s][ZONE1] = -(1 + config.compon[s].modulus.mRac);
      Frontiere[s][ZONE2] = config.compon[s].modulus.mRac;
      Frontiere[s][ZONE3] = config.compon[s].modulus.m - 1;
   }
   SetSeed (config.seed);
}

//===========================================================================

void Zone::init (const Component & comp, int j, int i)
{
   MScal h, Eb, Ec;

   // Calcul des bornes pour a[j][i] ou q (ZONE 1 ou 3) dans chaque zone
   if (comp.searchMethod != EXHAUST) {
      if (i == comp.k)
         h = comp.Hk;
      else
         h = comp.H;         // Taille des regions pour a[j][i]
   }
   Eb = comp.b[i];
   Ec = comp.c[i];

   if (comp.implemCond != APP_FACT) {
      // Pas de contrainte d'implantation; une seule grande zone
      inf = Eb;
      sup = Ec;                // Les bornes de la zone
      if (comp.searchMethod != EXHAUST) {
         // Ajuster taille des régions qui seront tirées au hasard
         frac = 1.0;
         T1 = sup - inf;
         smallF = T1 < h;
         if (!smallF)
            supMsH = sup - h + 1;
      }
      nextZone = 0;

   } else {  // La contrainte AppFact est appliquee; Max de 3 zones
      Zone *Z = 0, *pZ = 0;
      MScal Tot;
      Tot = 0;
      bool isFirst = true;
      int k = ZONE1;
      // Trouver la premiere zone k intersectant [Eb, Ec]
      while (k < NZONES - 1 && Frontiere[j][k] < Eb)
         ++k;

      // Initialisation de la zone k
      while (Eb <= Ec) {
         if (isFirst) {
            isFirst = false;
            Z = this;
         } else {
            // pZ precedes Z in the list of zones
            pZ->nextZone = Z = new Zone ();
            // setFrontieres (comp.modulus);
            if (comp.searchMethod != EXHAUST)
               ; // comp.ApproxTotGen = true;
         }
         pZ = Z;
         Z->No = (ZoneType) k;
         Z->calcInfBound (Eb, Ec, comp, k, Z->inf); // Bornes de la zone
         Z->calcSupBound (Eb, Ec, comp, k, Z->sup);
         T1 = Z->sup - Z->inf;        // Taille de la zone, -1
         Tot += T1;
         ++Tot;
         if (comp.searchMethod != EXHAUST) {
            // Ajuster taille des régions qui seront tirées au hasard
            Z->smallF = T1 < h;
            if (Z->smallF) {
               cout << "** WARNING ** Multiplier " << i <<
               " zone " << k << " : sup - inf + 1 < H(k). " << endl;
            } else {
               Z->supMsH = Z->sup - h + 1;
            }
         }
         Eb = Frontiere[j][k] + 1; // passer a la zone suivante
         ++k;
      }                        // while
      Z->nextZone = 0;

      // Calculer la proportion frac que represente chaque zone
      Z = this;
      double LRE, TotLR;
      while (Z != 0) {
         T1 = Z->sup - Z->inf + 1;
         conv (LRE, T1);
         conv (TotLR, Tot);
         Z->frac = LRE / TotLR;
         Z = Z->nextZone;
      }
   }
}


//===========================================================================

void Zone::calcInfBound (const MScal & b, const MScal & c,
                         const Component & comp, int No, MScal & inf)
{
   /* Calcule la borne inferieure de l'intersection d'une zone avec [b..c].
      Si No = 1 ou 3, la borne retournee sera (approx) q t.q. [m/q] = c
      Si No = 2, retourne b. A l'appel, b est toujours >= borne inf. de
      la zone. */

   switch (No) {
   case ZONE1:
      if (c < comp.modulus.mRacNeg) {
         T1 = comp.modulus.m - 1;
         T2 = c + 1;
         Divide (T1, T2, inf, T1);
         if (!IsZero (T1))
            inf += 1;
         if (inf < (comp.modulus.mRacNeg))
            inf = comp.modulus.mRacNeg;
      } else
         inf = comp.modulus.mRacNeg;
      break;

   case ZONE2:                      // -sqrt(m)-1 < b <= sqrt(m);
      inf = b;
      break;

   case ZONE3:
      /* b > sqrt(m) On retourne le a minimal tel que b <= floor(m/a) <= c.
         Note: floor(m/Borne) <= c, m/Borne < c+1, Borne > m/(c+1). On fait: 
         Borne = (m DIV (c+1)) + 1. */
      T1 = c + 1;
      Quotient (comp.modulus.m, T1, T2);
      inf = T2 + 1;
      if (inf == 1)               // si inf = 1, alors a = m
         inf = 2;
      break;

   default:
      cout << "Zone::calcInfBound:   impossible case" << endl;
      assert (0);
   }
}


//===========================================================================

void Zone::calcSupBound (const MScal & b, const MScal & c,
                         const Component & comp, int No, MScal & sup)
{
   switch (No) {
   case ZONE1:                      // b <= -sqrt(m)-1
      Quotient (comp.modulus.m, b, sup);
      if (sup == -1)               // si sup = -1, alors a = -m
         sup = -2;
      break;

   case ZONE2:                      // -sqrt(m)-1 < b <= sqrt(m)
      if (comp.modulus.mRac < c)
         sup = comp.modulus.mRac;
      else
         sup = c;
      break;

   case ZONE3:                      /* b > sqrt(m) On retourne le a maximal tel
                                                 que b <= floor(m/a) <= c. */
      Quotient (comp.modulus.m, b, sup);
      if (comp.modulus.mRac < sup)
         sup = comp.modulus.mRac;
      break;

   default:
      cout << "Zone::calcSupBound:   impossible case" << endl;
      assert (0);
   }
}


//===========================================================================

std::string Zone::toString()
{
   std::ostringstream sortie;
   sortie << "No = " << No << endl
   << "inf = " << inf << endl
   << "sup = " << sup << endl
   << "frac = " << frac << endl
   << "supMsH = " << supMsH << endl
   << "smallF = " << smallF << endl
   << "DivQ = { " << DivQ[0] << ", " << DivQ[1] << ", " <<
   DivQ[2] << " }" << endl
   << "Frontiere = { " << Frontiere[0][0] << ", " << Frontiere[0][1] << ", " <<
   Frontiere[0][2] << " }" << endl
   << nextZone << endl << endl;

   return sortie.str ();
}


//===========================================================================
/*
void Zone::toRegion (Zone *Z)
{
   inf = Z->getInf();
   sup = Z->getSup();
   No = Z->getNo ();
   frac = Z->getFrac();
   smallF = Z->smallF;
   supMsH = Z->getSupMsH();
   cond = Z->cond;
   nextZone = Z->nextZone;
   Frontiere[0] = Z->Frontiere[0];
   Frontiere[1] = Z->Frontiere[1];
   Frontiere[2] = Z->Frontiere[2];
}
*/

//===========================================================================

void Zone::chooseBoundaries (const Component & comp, Zone *zone)
{
   /* Choisir une region au hasard et en initialiser les bornes.
      On fouillera ensuite completement cette region.
      Utilise avec la methode de recherche RANDOM seulement. */

   Zone *Z;
   double p, u;
   Zone *R;
   long h;

   for (int i = 1; i <= comp.k; ++i) {
      Z = i + zone;
      if (Z->nextZone != 0) {
         // Choisir d'abord l'une des zones selon les probabilites Z->frac
         u = RandU01();
         p = Z->getFrac();
         while ((u > p) && (Z != 0)) {
            Z = Z->nextZone;
            p += Z->getFrac();
         }
      }

      // Z pointe a la zone choisie
      R = i + this;
      if (Z->smallF) {
         // On prends toute la zone
         R->getInf() = Z->getInf();
         R->getSup() = Z->getSup();
      } else {
         // On prends un intervalle de taille h, au hasard, dans la zone
         if (i == comp.k)
            h = comp.Hk - 1;
         else
            h = comp.H - 1;
         T1 = Z->getSupMsH() - Z->getInf() + 1;
         conv (p, T1);
         p *= RandU01();
         conv (T1, p);
         T1 = Z->getInf() + T1;
         R->getInf() = T1;
         R->getSup() = T1 + h;
      }
      R->No = Z->getNo();
   }
}


//===========================================================================

}
