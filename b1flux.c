/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : b1flux.c                                       */
/*                                                                           */
/* Created:       2013/04/14 (JLe)                                           */
/* Last modified: 2015/03/09 (JLe)                                           */
/* Version:       2.1.23                                                     */
/*                                                                           */
/* Description: Calculates B1 leakage-corrected flux                         */
/*                                                                           */
/* Comments:  - NumericGauss() vie about kaiken CPU-ajan --> optimointi      */
/*                                                                           */
/*            - 22.8.2014 / 2.1.22: NumericGauss() palauttaa -1 jos matriisi */
/*              singulaarinen. Tämä keskeyttää nyt iteroinnin. Jos noissa    */
/*              CRAM:istä lainatuissa aliohjelmissa varataan jotain          */
/*              muistia niin se jää siinä tapauksessa vapauttamatta. Koko    */
/*              rutiini pitää jossain vaiheessa kirjoittaa uusiksi.          */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "B1Flux:"

/* Setting maximum iterations to zero should reproduce infinite */
/* spectrum results. */

#define MAX_ITER 500
#define B2_LIM 1E-6

/*****************************************************************************/

long B1Flux(long gcu, long nmg, const double *flx0, const double *tot0, 
            const double *abs, const double *nsf, const double *sp0, 
            const double *sp1, const double *chit, double *flx1, double *D1)
{
  long n, m, i, iter, nnz, ptr;
  double *x, *alpha, sum1, sum2, L2, B2, val, kb1, d1, B, keff, err, J;
  double div, *tot;
  complex theta0, *ei, *di, **D, *chi, *flx,*LUval, *chival;
  struct ccsMatrix *Dinv, *A; 

  /***************************************************************************/

  /***** Calculate initial guess for keff and B2 *****************************/

  /* Reset convergence flag (to avoid including zero results in average) */

  WDB[gcu + GCU_B1_CONV] = (double)NO;

  /* Calculate and loss terms */

  sum1 = 0.0;
  sum2 = 0.0;
  
  for (n = 0; n < nmg; n++)
    {
      sum1 = sum1 + flx0[n]*nsf[n];
      sum2 = sum2 + flx0[n]*abs[n];
    }

  /* Check zero source rate (possible if B1 calculation is on and */
  /* universe has no fissile material) */

  if (sum1 == 0.0)
    {
      /* Set convergence flag (tämä mahdollistaa sen että vakiot  */
      /* kondensoidaan INF-spektrillä tapauksissa joissa lasku on */
      /* tuomittu epäonnistumaan). */
      
      WDB[gcu + GCU_B1_CONV] = (double)YES;
      
      /* Print warning */

      Note(0, "B1 calculation skipped due to zero source rate");
           
      /* Exit subroutine */

      return NO;
    }

  /* Check values */
  
  CheckValue(FUNCTION_NAME, "sum1", "", sum1, ZERO, INFTY);
  CheckValue(FUNCTION_NAME, "sum2", "", sum2, ZERO, INFTY);

  /* Calculate k-inf */
  
  keff = sum1/sum2;

  /* Skip B1 calculation if k-eff is too low (30.6.2014 / 2.1.22) */

  if (keff < 0.5)
    {
      /* Set convergence flag */
      
      WDB[gcu + GCU_B1_CONV] = (double)YES;
      
      /* Print warning */

      Note(0, "B1 calculation skipped due to low k-inf (%1.5f)", keff);
           
      /* Exit subroutine */

      return NO;
    }

  /* Add to counter */

  WDB[DATA_B1_REPEATED] = RDB[DATA_B1_REPEATED] + 1.0;

  /* Calculate transport cross section */

  sum1 = 0.0;
  sum2 = 0.0;
  
  for (n = 0; n < nmg; n++)
    {
      /* Average scattering cosine */
      
      val = 0.0;
      div = 0.0;
      
      for (m = 0; m < nmg; m++)
        val = val + flx0[n]*sp1[m*nmg + n];
      
      for (m = 0; m < nmg; m++)
        div = div + flx0[n]*sp0[m*nmg + n];
      
      /* Check */
      
      if (div != 0.0)
        val = val/div;
      else 
        val = 0.5;
      
      CheckValue(FUNCTION_NAME, "val", "", val, -1.0, 1.0);
      
      /* Add to transport cross section */
      
      if ((val = tot0[n] - val*(tot0[n] - abs[n])) > 0.0)
        sum1 = sum1 + flx0[n]/val;
      
      /* Add to flux */
      
      sum2 = sum2 + flx0[n];
    }
  
  /* Diffusion coefficient */
  
  if (sum2 > 0.0)
    d1 = sum1/sum2/3.0;
  else
    d1 = 1E-2;
  
  /* Calculate removal cross section (use absorption only) */
  
  sum1 = 0.0;
  sum2 = 0.0;
  
  for (n = 0; n < nmg; n++)
    {
      sum1 = sum1 + flx0[n]*abs[n];
      sum2 = sum2 + flx0[n];
    }
  
  /* Check values */
  
  CheckValue(FUNCTION_NAME, "sum1", "", sum1, ZERO, INFTY);
  CheckValue(FUNCTION_NAME, "sum2", "", sum2, ZERO, INFTY);
  
  /* Migration area */
  
  L2 = d1/(sum1/sum2);
  CheckValue(FUNCTION_NAME, "L2", "", L2, ZERO, INFTY);
  
  /* Buckling */
  
  B2 = (keff - 1.0)/L2;
  
  /*
  fprintf(out, "B2 = %11.5E, k = %1.5f \n", B2, keff);  
  */

  /***************************************************************************/

  /***** Allocate memory for temporary data **********************************/

  /* Allocate memory for temporary variables */

  x = (double *)Mem(MEM_ALLOC, nmg, sizeof(double));
  alpha = (double *)Mem(MEM_ALLOC, nmg, sizeof(double));

  /* Allocate memory and copy total xs to temporary array */

  tot = (double *)Mem(MEM_ALLOC, nmg, sizeof(double));
  memcpy(tot, tot0, nmg*sizeof(double));
  
  /* Ei tartte kopioida koko LU-matriisia vaan ainoastaan values joka kerta */

  LUval = (complex *)Mem(MEM_ALLOC, nmg*nmg, sizeof(complex));

  /* Allocate memory for matrix etc. */

  Dinv = ccsMatrixNew(nmg, nmg, nmg*nmg);
  A = ccsMatrixNew(nmg, nmg, nmg*nmg);

  ei = (complex *)Mem(MEM_ALLOC, nmg, sizeof(complex));   
  di = (complex *)Mem(MEM_ALLOC, nmg, sizeof(complex));   
  chi = (complex *)Mem(MEM_ALLOC, nmg, sizeof(complex));  
  chival = (complex *)Mem(MEM_ALLOC, nmg, sizeof(complex));  
  flx = (complex *)Mem(MEM_ALLOC, nmg, sizeof(complex));   

  /* Matriisi D kannattaa tallentaa suoraan taulukkona */

  D = (complex **)Mem(MEM_ALLOC, nmg, sizeof(complex *)); 

  for(n = 0; n < nmg; n++)
    D[n] = (complex *)Mem(MEM_ALLOC, nmg, sizeof(complex));

  /***************************************************************************/

  /***** Main loop ***********************************************************/

  /* Fission spectrum */
  
  for (n = 0; n < nmg; n++)
    {
      chi[n].re = chit[n];
      chi[n].im = 0.0;
    }
  
  /* Avoid compiler warning */

  kb1 = keff;
  
  /* Reset temporary variables */
  
  memset(x, 0.0, nmg*sizeof(double));
  memset(alpha, 0.0, nmg*sizeof(double));
  
  /* Reset error */
  
  err = 1.0;
  
  /* Iteration loop */
  
  /* Kaikki varsinainen laskenta tapahtuu tän silmukan sisällä.    */
  /* MAX_ITER on nyt asetetettu nollaksi, eli silmukkaa ei käydä   */
  /* ollenkaan läpi. Tuloksena saadut ryhmävakiot ovat tällöin     */
  /* NSF:ää lukuunottamatta numeroarvoltaan samat kuin infinite-   */
  /* spectrum -laskussa, mitä voi hyvin käyttää myöhemmin tehtävän */
  /* uudelleenhomogenisoinnin tarkastamiseen. Toi silmukan sisällä */
  /* tehtävä laskenta on vielä vähän mitä sattuu. Pääsin alkuun,   */
  /* mutta en saanut kaikista osista tolkkua. */

  for (iter = 0; iter < MAX_ITER; iter++)
    {
      /***********************************************************************/

      /***** Set temporary variables *****************************************/

      /* Check if too close to criticality */
          
      if (fabs(B2) < B2_LIM)
        break; 
        
      /* Square root of buckling */

      if (B2 >= 0.0)
        B = sqrt(B2);
      else
        B = sqrt(-B2);

      /* Calculate temporary variable x (toi voi mennä INF:ksi!) */

      for (n = 0; n < nmg; n++)
        {              
          x[n] = B/tot[n];

          /* NOTE: Tää voi mennä huonon statistiikan takia yli 1,      */
          /* mikä aiheuttaa NAN:in tuolla myöhemmin jos B2 < 0.        */
          /* Toi yläraja muutettiin 1 -> 0.999999 (1.12.2014 / 2.1.23) */
          
          if (x[n] > 0.999999)
            {
              tot[n] = B;
              x[n] = 0.999999;
            }
          
          CheckValue(FUNCTION_NAME, "x[n]", "", x[n], 0.0, 1.0);
        }       

      /* Tarkistetaan jos luuppi katkaistu, asetetaan iteraatioiden */
      /* määrä maksimiin --> spektriä ei päivitetä, ja katkaistaan  */
      /* iteraatio (JLE 2.1.4 / 3.4.2012) */
      
      if (n < nmg)
        {
          iter = MAX_ITER;
          break;
        }
        
      /* Calculate alpha  */
          
      if (B2 > 0.0)
        {            
          for (n = 0; n < nmg; n++)
            {
              alpha[n] = x[n]*x[n];
              alpha[n] = alpha[n]*atan(x[n])/(x[n] - atan(x[n])); 
              CheckValue(FUNCTION_NAME, "alpha[n] (1)", "", alpha[n], 
                         -INFTY, INFTY);
            }
        }
      else if(B2 < 0.0)
        {
          for (n = 0; n < nmg; n++)
            {
              alpha[n] = x[n]*x[n];
              alpha[n] = alpha[n]*log((1.0 + x[n])/(1.0 - x[n]))/ 
                (log((1.0 + x[n])/(1.0 - x[n])) - 2.0*x[n]);
                  
              CheckValue(FUNCTION_NAME, "alpha[n] (2)", "", alpha[n], 
                         -INFTY, INFTY);
            }
        }
      else
        {
          for (n = 0; n < nmg; n++)
            alpha[n] = 3.0; 
        }
      
      /***********************************************************************/

      /***** Form matrix Dinv ************************************************/

      /* Reset size and pointer to first column */

      nnz = 0;
      Dinv->colptr[0] = 0 ;

      /* Loop over columns */
      
      for (m = 0; m < nmg; m++)
        {
          /* Loop over rows */
          
          for (n = 0; n < nmg; n++)
            {
              /* Index to scattering matrixes */

              i = m*nmg + n;

              /* Add value */
              
              CheckValue(FUNCTION_NAME, "sp1", "", sp1[i], -INFTY, INFTY);
              Dinv->values[nnz].re = -3.0*sp1[i];
              
              /* Add diagonal */
              
              if (n == m)
                {
                  CheckValue(FUNCTION_NAME, "alpha", "", alpha[n], -INFTY, 
                             INFTY);
                  CheckValue(FUNCTION_NAME, "tot", "", tot[m], -INFTY, INFTY);
                  Dinv->values[nnz].re = Dinv->values[nnz].re +
                    alpha[n]*tot[m];
                }
              
              /* Reset imaginary part */
              
              Dinv->values[nnz].im = 0.0;
              
              CheckValue(FUNCTION_NAME, "Dinv->values[nnz].re", "", 
                         Dinv->values[nnz].re, -INFTY, INFTY);
              CheckValue(FUNCTION_NAME, "Dinv->values[nnz].im", "", 
                         Dinv->values[nnz].im, 0.0, 0.0);
              
              /* Put row index */
              
              Dinv->rowind[nnz] = n;
              
              /* Update size */
              
              nnz++;
            }        
          
          /* Put column index */
          
          Dinv->colptr[m + 1] = nnz;
        }

      /***********************************************************************/
          
      /***** Form matrix D ***************************************************/

#ifdef DEBUG
      
      for (i=0; i < Dinv->nnz; i++)
        {
          CheckValue(FUNCTION_NAME, "Dinv->values[i].re", "",  
                     Dinv->values[i].re, -INFTY, INFTY);
          CheckValue(FUNCTION_NAME, "Dinv->values[i].im", "", 
                     Dinv->values[i].im, 0.0, 0.0);
        }

#endif
      
      /* Ota values talteen */
      
      memcpy(LUval, Dinv->values, (Dinv->nnz) * sizeof(complex));

#ifdef DEBUG
      
      for (i=0; i < Dinv->nnz; i++)
        {
          CheckValue(FUNCTION_NAME, "LUval.re", "", LUval[i].re, -INFTY, INFTY);
          CheckValue(FUNCTION_NAME, "LUval.im", "", LUval[i].im, 0.0, 0.0);
        }

#endif

      /* Set theta to zero */
      
      theta0.re = 0.0; 
      theta0.im = 0.0; 

      FindRowIndexes(Dinv);

      /* Loop over energy groups */
      /*
      ResetTimer(TIMER_MISC);
      StartTimer(TIMER_MISC);
      */
      for (n = 0; n < nmg; n++)
        {
          /* Form ei */
          
          memset(ei, 0.0, nmg*sizeof(complex));
          ei[n].re = 1.0; 
          
          memcpy(Dinv->values, LUval, Dinv->nnz*sizeof(complex));
          
          /* D:n sarakkeille varattu tila valmiiksi, Gauss ei */
          /* varaa mitään ! */

          if ((NumericGauss(Dinv, ei, theta0, D[n])) < 0)
            {
              /* Free allocated memory */

              Mem(MEM_FREE, x);
              Mem(MEM_FREE, alpha);
              Mem(MEM_FREE, LUval);
              
              ccsMatrixFree(Dinv);
              ccsMatrixFree(A);
              
              Mem(MEM_FREE, ei); 
              Mem(MEM_FREE, di); 
              Mem(MEM_FREE, chi);
              Mem(MEM_FREE, chival);
              Mem(MEM_FREE, flx); 
              
              for (n = 0; n < nmg; n++)
                Mem(MEM_FREE, D[n]);
              
              Mem(MEM_FREE, D); 
              Mem(MEM_FREE, tot);
              
              /* Print warning */
              
              Note(0, "Error in B1 calculation");

              /* Exit subroutine */
              
              return NO;
            }
        }

      /*
      StopTimer(TIMER_MISC);
      printf("%f\n", TimerCPUVal(TIMER_MISC));
      */
      /***********************************************************************/
          
      /***** Form matrix A ***************************************************/

      /* Reset size and pointer to first column */
          
      nnz = 0;
      A->colptr[0] = 0 ;
      
      /* Loop over columns */
      
      for (m = 0; m < nmg; m++)
        {
          /* Loop over rows */
          
          for (n = 0; n < nmg; n++)
            {
              /* Index to scattering matrixes */

              i = m*nmg + n;

              /* Set value */
                  
              A->values[nnz].re = -sp0[i] + B2*D[m][n].re; /* D^T*/
                  
              /* Add diagonal */
              
              if (n == m)
                A->values[nnz].re = A->values[nnz].re + tot[n];
                  
              /* Reset imaginary part */
              
              A->values[nnz].im = 0.0;
              
              /* Put row index */
              
              A->rowind[nnz] = n;
              
              /* Update size */
              
              nnz++;
            }

          /* Put column index */
              
          A->colptr[m + 1] = nnz;
        }
          
      /* Set size */
      
      A->nnz = nnz;

      /***********************************************************************/
          
      /***** Solve equation **************************************************/
          
      /* Solve flux */
          
      FindRowIndexes(A);
      memcpy(chival, chi, nmg*sizeof(complex));

      /* HUOM! Myös chival gaussataan... */ 
    
      if (NumericGauss(A, chival, theta0, flx) < 0)
        {
          /* Free allocated memory */

          Mem(MEM_FREE, x);
          Mem(MEM_FREE, alpha);
          Mem(MEM_FREE, LUval);
          
          ccsMatrixFree(Dinv);
          ccsMatrixFree(A);
          
          Mem(MEM_FREE, ei); 
          Mem(MEM_FREE, di); 
          Mem(MEM_FREE, chi);
          Mem(MEM_FREE, chival);
          Mem(MEM_FREE, flx); 
          
          for (n = 0; n < nmg; n++)
            Mem(MEM_FREE, D[n]);
          
          Mem(MEM_FREE, D); 
          Mem(MEM_FREE, tot);
          
          /* Print warning */
          
          Note(0, "Error in B1 calculation");

          /* Exit subroutine */
          
          return NO;
        }
      
      /* B1 keff */
      
      kb1 = 0.0;
      for (n = 0; n < nmg; n++)
        kb1 = kb1 + nsf[n]*flx[n].re; 
      
      /* Interpolointi */
          
      val = B2/(1.0/kb1 - 1.0/keff); 
      B2 = B2 + val - val/kb1; 
      
      /*
      fprintf(out, "B2 = %11.5E, k = %1.5f \n", B2, kb1);  
      */

      /* Calculate error */
      
      err = kb1 - 1.0; 
      
      /* Check convergence */
      
      if (fabs(err) < RDB[DATA_B1_ERR_LIMIT])
        break; 

      /***********************************************************************/
    }

  /***************************************************************************/

  /***** Iteration completed *************************************************/

  /* Check convergence */
      
  if (iter < MAX_ITER)
    {
      /* Normalize critical spectrum */

      sum1 = 0.0;
      sum2 = 0.0;

      for (n = 0; n < nmg; n++)
        {
          sum1 = sum1 + flx0[n];
          sum2 = sum2 + flx[n].re;
        }

      /* Check */

      if (sum1 < ZERO)
        Die(FUNCTION_NAME, "Error in normalization");
      else if (sum2 < ZERO)
        {
          /* Free allocated memory */

          Mem(MEM_FREE, x);
          Mem(MEM_FREE, alpha);
          Mem(MEM_FREE, LUval);
          
          ccsMatrixFree(Dinv);
          ccsMatrixFree(A);
          
          Mem(MEM_FREE, ei); 
          Mem(MEM_FREE, di); 
          Mem(MEM_FREE, chi);
          Mem(MEM_FREE, chival);
          Mem(MEM_FREE, flx); 
          
          for (n = 0; n < nmg; n++)
            Mem(MEM_FREE, D[n]);
          
          Mem(MEM_FREE, D); 
          Mem(MEM_FREE, tot);
          
          /* Print warning */
              
          Note(0, "Error in B1 calculation");

          /* Exit subroutine */
          
          return NO;
        }

      /* Calculate normalized flux */

      for (n = 0; n < nmg; n++)
        flx1[n] = flx[n].re*sum1/sum2;

      /* Calculate B1 current and diffusion coefficient */

      for (n = 0; n < nmg; n++)
        {
          /* Current */

          J = 0.0;
          for (m = 0; m < nmg; m++)
            J = J + D[m][n].re*flx1[m];  
      
          /* Diffusion coefficient */
      
          if (flx1[n] > 0.0)
            D1[n] = J/flx1[n];
          else
            D1[n] = 0.0;
        }

      /* Store k-inf, B2 and error*/

      ptr = (long)RDB[gcu + GCU_B1_KEFF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddStat(kb1, ptr, 0);

      ptr = (long)RDB[gcu + GCU_B1_B2];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddStat(B2, ptr, 0);

      ptr = (long)RDB[gcu + GCU_B1_ERR];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddStat(fabs(err), ptr, 0);
      
      /* Set convergence flag */

      WDB[gcu + GCU_B1_CONV] = (double)YES;

      /* Add to counter */

      WDB[DATA_B1_CONVERGED] = RDB[DATA_B1_CONVERGED] + 1.0;
    }
  else
    {
      /* No convergence */
          
      Note(0, "No convergence in B1 iteration, final error = %1.5E.\n\n", 
           fabs(err));
    }
  
  /***************************************************************************/

  /*
  fprintf(out, "\nB1 iteration completed, final error %1.2E.\n\n", err);  
  */

  /* Free allocated memory */

  Mem(MEM_FREE, x);
  Mem(MEM_FREE, alpha);
  Mem(MEM_FREE, LUval);

  ccsMatrixFree(Dinv);
  ccsMatrixFree(A);

  Mem(MEM_FREE, ei); 
  Mem(MEM_FREE, di); 
  Mem(MEM_FREE, chi);
  Mem(MEM_FREE, chival);
  Mem(MEM_FREE, flx); 

  for (n = 0; n < nmg; n++)
    Mem(MEM_FREE, D[n]);

  Mem(MEM_FREE, D); 
  Mem(MEM_FREE, tot);

  /* Exit subroutine */

  if (iter < MAX_ITER)
    return YES;
  else
    return NO;
}

/*****************************************************************************/

