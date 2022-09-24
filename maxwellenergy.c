/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : maxwellenergy.c                                */
/*                                                                           */
/* Created:       2011/03/02 (JLe)                                           */
/* Last modified: 2011/10/28 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Samples energy from a Maxwellian distribution                */
/*                                                                           */
/* Comments: - From Serpent 1.1.0                                            */
/*           - Käytetään MCNP:n menetelmää. Lux & Koblinger antaa eri        */
/*             tuloksia kuin nopeuden komponenttien arpominen normaali-      */
/*             jakaumasta.                                                   */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "MaxwellEnergy:"

/*****************************************************************************/

double MaxwellEnergy(double kT, long id)
{
  double rand1, rand2, rand3, R;
  double rand4;

  /* Use the method in MCNP4c (manual p. 2-44) */

  do
    {
      rand1 = RandF(id);
      rand2 = RandF(id);
      R = rand1*rand1 + rand2*rand2;
    } 
  while(R > 1.0);
  
  rand3 = log(RandF(id));
  rand4 = log(RandF(id));

  return -kT*(rand1*rand1*rand3/R + rand4);
}

/*****************************************************************************/
