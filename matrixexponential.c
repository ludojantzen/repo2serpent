/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : matrixexponential.c                            */
/*                                                                           */
/* Created:       2011/05/03 (MPu)                                           */
/* Last modified: 2015/04/08 (JLe)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description:                                                              */
/*                                                                           */
/* Comments: Funktio palauttaa matriisieksponenttiratkaisun N=exp(At)N0,     */
/*           missä exp(At)N0 on laskettu CRAM-menetelmän kertaluvulla k      */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MatrixExponential:"

/*****************************************************************************/

/* A = palamamatriisi, N0 = nukliditiehydet alussa, t = aika-askel,   */
/* cram_k = CRAM-kertaluku, funktio palauttaa ratkaisun N = exp(At)N0 */

double *MatrixExponential(struct ccsMatrix *A, double *N0, double t)
{
  long i, j, k, n, m, cram_k;
 
  double *N;  /* CRAM-ratkaisu, joka palautetaan */

  complex z, alpha0, alpha[8], theta[8], *b, *N0c, *x, *y, *val; 

  struct ccsMatrix *LU;

  /***************************************************************************/
  /* Ota valitun kertaluvun mukaiset CRAM-kertoimet                          */

  /* k = 4, 6, 8, 10, 12, 14 tai 16 */
  /* Otin tosta muistinvarauksen pois ja muutin alphan ja thetan */
  /* kiinteän pituisiksi taulukoiksi 5.6.2012 / 2.1.6 (JLe) */

  cram_k = (long)RDB[DATA_BURN_CRAM_K];

  switch(cram_k)
    {
    case 4:
      {
        alpha0.re  = +8.647853972498786800e-005;
        alpha0.im  = +0.000000000000000000e+000;  
        
        alpha[ 0].re = -6.168352255495653700e-002; 
        alpha[ 1].re = +7.339241923416284200e-002; 
        alpha[ 0].im = -1.905059455979927400e+000; 
        alpha[ 1].im = +4.500049158538028100e-001; 
        
        theta[ 0].re = +1.548400570539430600e+000; 
        theta[ 1].re = -3.678383143998292200e-001; 
        theta[ 0].im = +1.191825853927608500e+000; 
        theta[ 1].im = +3.658133272063312600e+000; 
        
        break; 
      }
    case 6: 
      {
        alpha0.re  = +1.008500923882493800e-006; 
        alpha0.im  = +0.000000000000000000e+000; 
        
        alpha[ 0].re = -5.790130040315515400e-001; 
        alpha[ 1].re = +6.630068705292500000e-001; 
        alpha[ 2].re = -8.358161715649581300e-002; 
        alpha[ 0].im = -4.286888564579144600e+000; 
        alpha[ 1].im = +1.451412919901465800e+000; 
        alpha[ 2].im = -1.064292607477132400e-001; 
        
        theta[ 0].re = +2.400602938932869400e+000; 
        theta[ 1].re = +1.158552571719206200e+000; 
        theta[ 2].re = -1.781988275920925300e+000; 
        theta[ 0].im = +1.193129308401788400e+000; 
        theta[ 1].im = +3.614772600818597000e+000; 
        theta[ 2].im = +6.196512467344857300e+000; 
        
        break; 
      }
    case 8:
      {
        alpha0.re  = +1.172260916115774900e-008; 
        alpha0.im  = +0.000000000000000000e+000; 
        
        alpha[ 0].re = -1.831771710726351800e+000; 
        alpha[ 1].re = +2.436240733060475600e+000; 
        alpha[ 2].re = -6.325880537788576300e-001; 
        alpha[ 3].re = +2.812975716585542500e-002; 
        alpha[ 0].im = -9.525608129487752200e+000; 
        alpha[ 1].im = +3.716755640605497600e+000; 
        alpha[ 2].im = -4.439231026351748700e-001; 
        alpha[ 3].im = +1.157738456447322300e-002; 
        
        theta[ 0].re = +3.220945245025599800e+000; 
        theta[ 1].re = +2.292249147806407500e+000; 
        theta[ 2].re = +2.694909873896886300e-001; 
        theta[ 3].re = -3.408539501369696600e+000; 
        theta[ 0].im = +1.193619605400713100e+000; 
        theta[ 1].im = +3.600771496019942400e+000; 
        theta[ 2].im = +6.082032592594966800e+000; 
        theta[ 3].im = +8.773034564209789100e+000;
        
        break;  
      }
    case 10:
      {
        alpha0.re  = +1.361066814808964400e-010; 
        alpha0.im  = +0.000000000000000000e+000; 
        
        alpha[ 0].re = -4.818382049919602400e+000; 
        alpha[ 1].re = +7.117165194845109800e+000; 
        alpha[ 2].re = -2.565584990860251000e+000; 
        alpha[ 3].re = +2.725869847632842500e-001; 
        alpha[ 4].re = -5.784903977272644800e-003; 
        alpha[ 0].im = -2.105459756329741000e+001; 
        alpha[ 1].im = +8.819533309589196800e+000; 
        alpha[ 2].im = -1.216385735471696700e+000; 
        alpha[ 3].im = +1.421172898473725200e-002; 
        alpha[ 4].im = +6.858506635486216500e-004; 
        
        theta[ 0].re = +4.027732482085530000e+000; 
        theta[ 1].re = +3.283752898170031000e+000; 
        theta[ 2].re = +1.715406031280785900e+000; 
        theta[ 3].re = -8.944046850939565400e-001; 
        theta[ 4].re = -5.161191252114784700e+000; 
        theta[ 0].im = +1.193856067339377700e+000; 
        theta[ 1].im = +3.594386774977728700e+000; 
        theta[ 2].im = +6.038934929965273900e+000; 
        theta[ 3].im = +8.582756905792797000e+000; 
        theta[ 4].im = +1.137515626304168000e+001; 
        
        break; 
      }
    case 12:
      {
        alpha0.re  = +1.591948795010012000e-012; 
        alpha0.im  = +0.000000000000000000e+000; 
        
        alpha[ 0].re = -1.179938501315198400e+001; 
        alpha[ 1].re = +1.878598503026608300e+001; 
        alpha[ 2].re = -8.238258735818444500e+000; 
        alpha[ 3].re = +1.319411773677646700e+000; 
        alpha[ 4].re = -6.857148341649022000e-002; 
        alpha[ 5].re = +8.184330422804526100e-004; 
        alpha[ 0].im = -4.641164503553920200e+001; 
        alpha[ 1].im = +2.023728672471359900e+001; 
        alpha[ 2].im = -2.796189543077298300e+000; 
        alpha[ 3].im = -1.835242229388347200e-001; 
        alpha[ 4].im = +3.841913592313700400e-002; 
        alpha[ 5].im = -5.813543637533003400e-004; 
        
        theta[ 0].re = +4.827493721838796300e+000; 
        theta[ 1].re = +4.206124482219006900e+000; 
        theta[ 2].re = +2.917868831275860100e+000; 
        theta[ 3].re = +8.517073748131417400e-001; 
        theta[ 4].re = -2.235968024316078200e+000; 
        theta[ 5].re = -6.998687843699482200e+000; 
        theta[ 0].im = +1.193987934738107300e+000; 
        theta[ 1].im = +3.590920581768813700e+000; 
        theta[ 2].im = +6.017345604827832000e+000; 
        theta[ 3].im = +8.503832330996919700e+000; 
        theta[ 4].im = +1.110929552413979600e+001; 
        theta[ 5].im = +1.399591568135529400e+001; 
        
        break; 
      }
    case 14: 
      {
        alpha0.re  = +1.832174378254041200e-014; 
        alpha0.im  = +0.000000000000000000e+000; 
        
        alpha[ 0].re = -7.154288063589066700e-005; 
        alpha[ 1].re = +9.439025310736169600e-003; 
        alpha[ 2].re = -3.763600387822697000e-001; 
        alpha[ 3].re = -2.349823209108270200e+001; 
        alpha[ 4].re = +4.693327448883129400e+001; 
        alpha[ 5].re = -2.787516194014564500e+001; 
        alpha[ 6].re = +4.807112098832508500e+000; 
        alpha[ 0].im = +1.436104334954130000e-004; 
        alpha[ 1].im = -1.718479195848301700e-002; 
        alpha[ 2].im = +3.351834702945010200e-001; 
        alpha[ 3].im = -5.808359129714207600e+000; 
        alpha[ 4].im = +4.564364976882776400e+001; 
        alpha[ 5].im = -1.021473399905645200e+002; 
        alpha[ 6].im = -1.320979383742872500e+000; 
        
        theta[ 0].re = -8.897773186468889200e+000; 
        theta[ 1].re = -3.703275049423448100e+000; 
        theta[ 2].re = -2.087586382501301400e-001; 
        theta[ 3].re = +3.993369710578568400e+000; 
        theta[ 4].re = +5.089345060580624700e+000; 
        theta[ 5].re = +5.623142572745977400e+000; 
        theta[ 6].re = +2.269783829231112500e+000; 
        theta[ 0].im = +1.663098261990208500e+001; 
        theta[ 1].im = +1.365637187148326800e+001; 
        theta[ 2].im = +1.099126056190126000e+001; 
        theta[ 3].im = +6.004831642235037000e+000; 
        theta[ 4].im = +3.588824029027006400e+000; 
        theta[ 5].im = +1.194069046343966900e+000; 
        theta[ 6].im = +8.461737973040222000e+000; 
        
        break; 
      }
    case 16:
    default:
      {
        cram_k = 16;

        alpha0.re  = +2.124853710495223600e-016; 
        alpha0.im  = +0.000000000000000000e+000; 
        
        alpha[ 0].re = -5.090152186522492000e-007; 
        alpha[ 1].re = +2.115174218246603100e-004; 
        alpha[ 2].re = +1.133977517848393000e+002; 
        alpha[ 3].re = +1.505958527002346700e+001; 
        alpha[ 4].re = -6.450087802553964900e+001; 
        alpha[ 5].re = -1.479300711355799900e+000; 
        alpha[ 6].re = -6.251839246320791700e+001; 
        alpha[ 7].re = +4.102313683541002100e-002; 
        alpha[ 0].im = -2.422001765285228900e-005; 
        alpha[ 1].im = +4.389296964738067200e-003; 
        alpha[ 2].im = +1.019472170421585700e+002; 
        alpha[ 3].im = -5.751405277642182200e+000; 
        alpha[ 4].im = -2.245944076265209600e+002; 
        alpha[ 5].im = +1.768658832378293700e+000; 
        alpha[ 6].im = -1.119039109428322800e+001; 
        alpha[ 7].im = -1.574346617345546700e-001; 
        
        theta[ 0].re = -1.084391707869698800e+001; 
        theta[ 1].re = -5.264971343442646800e+000; 
        theta[ 2].re = +5.948152268951177200e+000; 
        theta[ 3].re = +3.509103608414918100e+000; 
        theta[ 4].re = +6.416177699099434600e+000; 
        theta[ 5].re = +1.419375897185666000e+000; 
        theta[ 6].re = +4.993174737717996700e+000; 
        theta[ 7].re = -1.413928462488886200e+000; 
        theta[ 0].im = +1.927744616718165100e+001; 
        theta[ 1].im = +1.622022147316792800e+001; 
        theta[ 2].im = +3.587457362018322100e+000; 
        theta[ 3].im = +8.436198985884374200e+000; 
        theta[ 4].im = +1.194122393370138600e+000; 
        theta[ 5].im = +1.092536348449672300e+001; 
        theta[ 6].im = +5.996881713603942300e+000; 
        theta[ 7].im = +1.349772569889274500e+001; 
                
        break; 
      }
    }
  
  /**************************************************************************/

  /* Aseta muuttujat kenttiin */

  n   = A->n; 
  m   = A->m; 
  
  if (m != n)
    Die(FUNCTION_NAME, "Burnup matrix not square\n"); 
  
  /* Varaa tilaa uusille muuttujille */
  /* Lasketaan ratkaisu ensin tähän */

  y =  (complex *)Mem(MEM_ALLOC, n, sizeof(complex));

  /* Reaalinen ratkaisu joka palautetaan */


  N = (double  *)Mem(MEM_ALLOC, n, sizeof(double));

  /*
  N = WorkArray(DATA_PTR_WORK_COMP2, PRIVA_ARRAY, n, OMP_THREAD_NUM);
  */

  /* b = alpha[k] * N0 */

  b = (complex *)Mem(MEM_ALLOC, n, sizeof(complex));

  /* (A - theta[k]*I) x = b */

  x = (complex *)Mem(MEM_ALLOC, n, sizeof(complex));

  /* N0 kompleksiformaatissa  */
  
  N0c = (complex *)Mem(MEM_ALLOC, n, sizeof(complex));

/*****************************************************************************/

  /* Tarkastetaan matriisin A alkiot */

#ifdef DEBUG

  for (i = 0; i < A->nnz; i++)
    {
      CheckValue(FUNCTION_NAME, "A re (1)", "", A->values[i].re, -INFTY, 
                 INFTY);
      CheckValue(FUNCTION_NAME, "A im (1)", "", A->values[i].im, 0.0, 0.0);
    }

#endif

  /* Lasketaan ensin symbolinen LU-hajotelma (tämä tarvitsee tehdä */
  /* vain kerran) */

  LU = SymbolicLU(A);
  /* Hae rivi-informaatio numeerista eliminaatiota varten */
  FindRowIndexes(LU);

  /* Tarkastetaan matriisin LU alkiot */
  
#ifdef DEBUG

  for (i = 0; i < LU->nnz; i++)
    {
      CheckValue(FUNCTION_NAME, "LU re (1)", "", LU->values[i].re, -INFTY, 
                 INFTY);
      CheckValue(FUNCTION_NAME, "LU im (1)", "", LU->values[i].im, 0.0, 0.0);
    }
  
#endif



  /* ccsMatrixPrint(LU); */
  
  /* Varataan tila toiselle samankokoiselle matriisille */
  
  /* Jaa, onkohan tässä turhaa kopioida koko matriisi, eikö riittäisi ottaa talteen values??*/

  /*LU_copy = ccsMatrixNew(n,n,LU->nnz);*/

  /* Tarkistetaan aika-askel (JLe) */

  if ((t < ZERO) || (t > INFTY))
    Die(FUNCTION_NAME, "t = %E\n", t);

  /* Alusta muuttujat  */

  for (i=0 ; i < LU->nnz; i++)
    LU->values[i].re = LU->values[i].re*t; /* Lisää aika t matriisiin */
  
  for (i=0; i < n; i++)
    N0c[i].re = N0[i]; 

  /* Ota matriisin numeeriset arvot talteen */

  val = (complex *)Mem(MEM_ALLOC, LU->nnz, sizeof(complex));
  memcpy(val, LU->values, (LU->nnz) * sizeof(complex)); 
    
  /***************************************************************************/

  /* Lasketaan CRAM-approksimaatio astetta cram_k */
  
  /* vaatii cram_k / 2 matriisin kääntöä */
  
  for (k=0; k < cram_k/2; k++)
    { 
      /* Kopioi LU -> LU_copy */

      memcpy(LU->values, val, (LU->nnz) * sizeof(complex)); 

      /*    ccsMatrixCopy(LU, LU_copy); */
      /* tarvitseeko tässä oikeestaan kopioida muuta kuin values? */
      /* Eikö indeksit ja muut pysy samana koko ajan?? */


      /* b = alpha(k) * N0 */
    
      for (j=0; j < n; j++)
        b[j] = c_mul( alpha[k], N0c[j]); 
    
      /* ratkaise (A - theta[k]*I)x = b  */

      NumericGauss(LU, b, theta[k], x); 
       
      /* y = alpha0 * N0 + sum_i (A - theta[k]*I)^(-1) * b */
       
      for (i= 0; i < n; i++)
        y[i] = c_add(y[i], x[i]);               
    }

  /* Lisätään vielä alpha0 * N0 muuttujaan y */
  
  for (j=0; j < n; j++)
    {
      z = c_mul(alpha0, N0c[j]);     
      N[j] = (2.0 * y[j].re) + z.re;
    }
  
  /**************************************************************************/
  
  /* Vapauta muisti */

  Mem(MEM_FREE, y);
  Mem(MEM_FREE, b);
  Mem(MEM_FREE, N0c);
  Mem(MEM_FREE, x, 0);
  ccsMatrixFree(LU);
  
  Mem(MEM_FREE, val);

  /*  ccsMatrixFree(LU_copy); */

  /* palauta ratkaisu */

  return N; 
  
  /**************************************************************************/
}

/*****************************************************************************/
