/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : ttaloop.c                                      */
/*                                                                           */
/* Created:       2011/06/10 (AIs)                                           */
/* Last modified: 2011/10/28 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Ratkaisu TTA-ketjuille, joissa on toistoja.                  */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "TTALoop:"

#define MAX_TTA_CHAIN 500  /* ketjun maksimipituus */

/* Kaksi hajoamisvakiota lasketaan yhtä suuriksi,
   jos suhteellinen ero on alle tämän. */

#define LIM_EQ 1.0E-5                                               

/* Variaation suuruus. Toistetut hajoamisvakiot kerrotaan tällä, 
   jotta niistä saadaan erillisiä. */

#define VARIATION_MULTIPLIER 1.0001

/* Katkaise ketju jos P kasvaa enemmän kuin tällä tekijällä.
   Antaa numeeriselle heitolle tilaa. */

#define LIM_P_STAB 1.001

/*****************************************************************************/

void TTALoop (long iso, double t, long n, double *l, double B, double Pin,
              double adens, double Plimit, double *N, struct ccsMatrix *A, 
              long initial)
{
  long i,j, dgt;
  double lambda, lambdaReal, Pout, f, br;

  /* Tarkastetaan ketjun pituus */

  if (n > MAX_TTA_CHAIN - 1)
    Die(FUNCTION_NAME, "Max length of TTA chain exceeded.\n");

  /* Ketjun leikkaus jäljelle jäävän osuuden perusteella */
  
  if (Pin < Plimit)
    return;
  
  /* transmutaatiokerroin */
  
  lambda=0.0;
  for (j=A->colptr[iso] ; j < A->colptr[iso+1] ; j++)
    if(A->rowind[j]==iso)
      {
        lambda=-A->values[j].re;
        break;
      }
  
  lambdaReal=lambda; /* käytetään branching ratioiden saamiseen */

  /* Jos nuklidi on stabiili... */
  
  if (lambda == 0.0)
    {
      N[iso] = N[iso] + Pin*adens;
      return;
    }

  /* Estetään mahdolliset toistot varioimalla */

  for(i = 0; i < n; i++)
    if (fabs((lambda - l[i]) / l[i]) < LIM_EQ)
      {
        lambda=lambda*VARIATION_MULTIPLIER;
        i=-1; /* uudestaan, uudestaan */
      }

  /* Lisätään nuklidi */
  
  l[n]=lambda;
  n++;
  
  /* Ratkaistaan ketju */
  
  f = TTAChain(n,t,l, &Pout);

  Pout=Pout*B;

  /* Tarkastetaan arvot */

  if (f > 1.0 || f < 0.0)
    return; /* numeriikka petti */
  
  if(Pout < 0.0 || Pout > Pin*LIM_P_STAB)
    return; /* numeriikka petti */

  /* Lisätään osuus tiheyteen */
  
  N[iso] = N[iso] + B*f*adens;

  /* looppi rivien, eli "tytärten yli"*/

  for (j=A->colptr[iso] ; j < A->colptr[iso+1] ; j++)
    {
      if(A->rowind[j]==iso)
        {
          continue; /* diagonaalialkio */
        }
      else
        {
          /* tytär jonka indeksi on */
          dgt=A->rowind[j];
          /* ja haarautumissuhde*/
          br=A->values[j].re/lambdaReal;
          
          /* Askel eteenpäin */
	  
          TTALoop(dgt, t, n, l, B*br, Pout*br, adens, Plimit, N, A, initial);
        }
    }
}

/*****************************************************************************/
