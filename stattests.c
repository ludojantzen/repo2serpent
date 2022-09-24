/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : stattests.c                                    */
/*                                                                           */
/* Created:       2013/06/12 (JLe)                                           */
/* Last modified: 2017/02/23 (TKa)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Subroutine for statistical tests                             */
/*                                                                           */
/* Comments: - Pääohjelma luuppaa statistiikkamuuttujien yli ja lukee datan  */
/*             vektoriin x. Varsinaiset testit voi tehdä paikalliseen        */
/*             aliohjelmaan StatTests0().                                    */
/*           - Bins are in different order in *_stat.m than in *_det*.m.     */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "StatTests:"

/* Local function definitions */

void StatTests0(double *, long, char *, FILE *);
static void PrintStattests(FILE *, long, double, double, long);
static long ShapiroWilkTest(double *, long, double *, double *);
static long DAgostinoPearsonTest(double *, long, double *, double *);
static long SkewnessTest(double *, long, double *, double *);
static long KurtosisTest(double *, long, double *, double *);
static long DriftInMeanTest(double *, long, double *, double *,long);
static long tCDF(double, double, double *);
static long tCDFNormalApprox(double, double, double *);
static long Betareg(double, double, double, double, long, double *);
static double SampleMean(double *, long);
static int CompareElements(const double *, const double *);
static double Polyval(const double *, long, double);
static double NormCDF(double);
static long NormCDFInv(double, double *);

/*****************************************************************************/

void StatTests()
{
  long loc0, ptr, i, n0, n1, n2, dim, bins[3], N;
  double *x;
  char tmpstr[MAX_STR], name[MAX_STR];
  FILE *fp;

  /* Avoid compiler warning */
  bins[0] = bins[1] = bins[2] = 0;

  /* Check if tests are run */

  if ((long)RDB[DATA_RUN_STAT_TESTS] == NO)
    return;
  
  /* Check mpi task (tarvitaan että rinnakkaismoodissa useampi taski */
  /* ei ala kirjoittamaan samaan fileen) */

  if (mpiid > 0)
    return;

  /* Set output file name */

  sprintf(tmpstr, "%s_stat.m", GetText(DATA_PTR_INPUT_FNAME));

  /* Open file */

  if ((fp = fopen(tmpstr, "w")) == NULL)
    Die(FUNCTION_NAME, "Unable to open file for writing");

  /* Print */

  fprintf(outp, "Performing statistical tests on output variables...\n");
  
  /* Print header to output file */
  
  fprintf(fp, "%% z = test statistic\n");
  fprintf(fp, "%% p = p-value calculated from z\n");
  fprintf(fp, "%% P = pass value (1 if p >= alpha, 0 if p < alpha)\n");
  fprintf(fp, "%%\n");
  fprintf(fp, "%% Used significance levels:\n");
  fprintf(fp, "%% Shapiro-Wilk        : alpha = 0.05\n");
  fprintf(fp, "%% D'Agostino-Pearson  : alpha = 0.05\n");
  fprintf(fp, "%% Drift-in-Mean       : alpha = 0.01\n");
  fprintf(fp, "%%\n");
  fprintf(fp, "%%%48s%34s%34s\n", "Shapiro-Wilk test", "D'Agostino-Pearson test", 
          "Drift-in-Mean test");
  fprintf(fp, "%%%14s", "");
  for (i = 0; i < 3; i++)
    fprintf(fp, "%16s%14s%4s", "z", "p", "P");
  fprintf(fp, "\n");
  

  /* Allocate memory for array of values (muistinhallintaan käytetään */
  /* Serpentin omia funktioita) */

  N = (long)RDB[DATA_CYCLE_IDX] + 1; /* DATA_CYCLE_IDX starts from 0 */
  x = (double *)Mem(MEM_ALLOC, N, sizeof(double));

  /* Loop over statistics (Serpentin perustietorakenne on linkitetty lista, */
  /* johon pääsee käsiksi pointtereiden kautta, ja listan alkioiden välillä */
  /* liikutaan NextItem() ja PrevItem() -funktioiden avulla. Kaikki statis- */
  /* tiikka tallennetaan listaan jonka ensimmäiseen alkioon haetaan tässä   */
  /* pointteri) */

  loc0 = (long)RDB[DATA_PTR_SCORE0];  
  while(loc0 > VALID_PTR)
    {
      /* Check history pointer (Muuttujille ei automaattisesti varata tilaa */
      /* sukupolvikohtaisten arvojen tallennusta varten. Tämä pitää tehdä   */
      /* erikseen funktiolla AllocStatHistory() samassa yhteydessä kun      */
      /* muuttuja määritellään (katso esim. processstats.c -aliohjelma)).   */

      if ((long)RDB[loc0 + SCORE_PTR_HIS] > 0)
        {
          /* Get dimension (Rakenne, johon data tallennetaan, voi       */
          /* periaatteessa sisältää mielivaltaisen määrän dimensioita   */
          /* parametrien (taulukointi-indeksit) lukumäärästä riippuen.) */

          if ((dim = (long)RDB[loc0 + SCORE_DIM]) > 3)
            Die(FUNCTION_NAME, "Need more dimensions");

          /* Pointer to list of maximum values (CheckPointer() -funktio on  */
          /* debuggausta varten. Jos lähdekoodi on käännetty -DDEBUG option */
          /* kanssa, pointterin arvo tarkistetaan ja koodi antaa virhe-     */
          /* ilmoituksen jos jotain on pielessä. Näillä tarkistuksilla      */
          /* löydetään suurin osa pointterivirheistä ennen kuin ne          */
          /* aiheuttavat muistivuotoja ja muuta vakavampaa) */

          ptr = (long)RDB[loc0 + SCORE_PTR_NMAX];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Read bin sizes (parametrin arvojen lukumäärä) */

          for (i = 0; i < dim; i++)
            {
              /* Get number of bins */

              bins[i] = (long)RDB[ptr++];
            }

          /* Get variable name (nimi joka on tallennettu kun muuttuja on */
          /* määritelty NewStat() -funktiolla) */

          sprintf(name, "%s", GetText(loc0 + SCORE_PTR_NAME));

          /* Check number of dimensions (HisVal()-funktio jolla tuo data     */
          /* luetaan ottaa argumenttina muuttujan dimension verran indeksejä */
          /* eikä datan käsittelyyn ole oikein mitään systemaattista tapaa.  */
          /* Tässä on nyt eri luupit 1-, 2- ja 3-dimensioisille muuttujille. */

          fprintf(fp, "\n");
          fprintf(fp, "STATTESTS_%s = [\n", name);

          if (dim == 1)
            {
              /* 1D variables */

              for (n0 = 0; n0 < bins[0]; n0++)
                {
                  /* Loop over cycles and read data */

                  for (i = 0; i < N; i++)
                    x[i] = HisVal(loc0, i, n0);

                  /* Print */

                  fprintf(fp, " %4ld %4d %4d", n0 + 1, 1, 1);

                  StatTests0(x, N, name, fp);

                }
            }
          else if (dim == 2)
            {
              /* 2D variables */

              for (n0 = 0; n0 < bins[0]; n0++)
                for (n1 = 0; n1 < bins[1]; n1++)
                  {
                    /* Loop over cycles and read data */

                    for (i = 0; i < N; i++)
                      x[i] = HisVal(loc0, i, n0, n1);

                    /* Print */

                    fprintf(fp, " %4ld %4ld %4d", n0 + 1, n1 + 1, 1);

                    StatTests0(x, N, name, fp);

                  }
            }
          else if (dim == 3)
            {
              /* 3D variables */

              for (n0 = 0; n0 < bins[0]; n0++)
                for (n1 = 0; n1 < bins[1]; n1++)
                  for (n2 = 0; n2 < bins[2]; n2++)
                    {
                      /* Loop over cycles and read data */

                      for (i = 0; i < N; i++)
                        x[i] = HisVal(loc0, i, n0, n1, n2);

                      /* Print */

                      fprintf(fp, " %4ld %4ld %4ld", n0 + 1, n1 + 1, n2 + 1);

                      StatTests0(x, N, name, fp);

                    }
            }
          else
            Die(FUNCTION_NAME, "Error in dimension");

          fprintf(fp, "];\n");
        }

      /* Next variable */

      loc0 = NextItem(loc0);

    }

  /* Free memory */

  Mem(MEM_FREE, x);
  
  /* Close file */

  fclose(fp);
  
  /* Print */

  fprintf(outp, "OK.\n\n");

  /***************************************************************************/
}

/*****************************************************************************/


/*****************************************************************************/
void StatTests0(double *x, long N, char *name, FILE *fp)
{ 
  /* Performs the statistical tests and prints the output.
   */
  
  /***** Variables for the statistical tests *********************************/

  double Z_W_R = -1, p_W_R = -1, K2 = -1, p_K2 = -1;
  double Z_drift2[1] = { -1 };
  double p_drift2[1] = { -1 };
  long retsw, retdp, retdrift2;

  /* Set the significance levels */

  double alpha_sw = 0.05;          /* Shapiro-Wilk test */
  double alpha_dp = 0.05;          /* D'Agostino-Pearson test */
  double alpha_drift = 0.01;       /* Drift-in-Mean test */

  /***************************************************************************/


  /***** Variables for calculating the number of batches etc *****************/

  long i, j, imax, Nbatch, Ni, batch_size, Ninactive, Nactive;
  double scoresum;
  double *xbatch, *xactive;

  /* Separate treatment for gcu parameters which batch size is given by
   * RDB[DATA_MICRO_CALC_BATCH_SIZE] */
  char *paramspre_gcu[4] = {"B1_", "INF_", "FG_", "DF_"};   /* Prefixes */
  int is_gcu;  /* Flag for gcu */

  /* Parameters which batch size is one */
  char *params_batchsize_1[4] = {"MEAN_POP_SIZE", "MEAN_POP_WGT",
    "TRANSPORT_RUNTIME", "TRANSPORT_CPU_USAGE"};
  int is_batchsize_1;      /* Flag for batch_size = 1 */

  /***************************************************************************/


  /***** Find out the correct number of batches  *****************************/

  /* Number of inactive and active cycles */

  Ninactive = (long) RDB[DATA_CRIT_SKIP];
  Nactive = (long) RDB[DATA_CRIT_CYCLES];
  
  /* Check N */
  
  if (N != Ninactive + Nactive)
    Die(FUNCTION_NAME, "Number of data (%ld) is different than the sum of "
        "inactive (%ld) and active cycles (%ld)", N, Ninactive, Nactive);

  /* Test only active cycles */

  xactive = x + Ninactive;

  /* Check if the batch size is always one */

  is_batchsize_1 = 0;
  imax = (long)(sizeof(params_batchsize_1)/sizeof(params_batchsize_1[0]));
  for (i = 0; i < imax; i++) {
    if (!strcmp(name, params_batchsize_1[i])) {
      is_batchsize_1 = 1;
      break;
    }
  }

  /* Check if the parameter is a gcu parameter */

  is_gcu = 0;
  imax = (long)(sizeof(paramspre_gcu)/sizeof(paramspre_gcu[0]));
  if (!is_batchsize_1) {
    for (i = 0; i < imax; i++) {
      if (!strncmp(name, paramspre_gcu[i], strlen(paramspre_gcu[i]))) {
        is_gcu = 1;
        break;
      }
    }
  }

  /* Get the batch size */

  if (is_gcu) {
    /* Batch size is the gcu batch size */
    batch_size = (long)RDB[DATA_MICRO_CALC_BATCH_SIZE];
  }
  else if (is_batchsize_1) {
    /* Batch size is always one */
    batch_size = 1;
  }
  else {
    /* Batch size is the batch interval set by "pop" */
    batch_size = (long)RDB[DATA_BATCH_INTERVAL];
  }

  /* Number of batches (assuming N mod batch_size = 0) */

  Nbatch = Nactive/batch_size;

  /* Check that the number of batches is positive */

  if (Nbatch < 1) {
    Warn(FUNCTION_NAME, "Variable %s has zero or negative number of batches, statistical tests are not performed.", name);
    return;
  }

  /* Allocate memory */

  xbatch = (double *)Mem(MEM_ALLOC, Nbatch, sizeof(double));


  /* Copy batch values to xbatch */

  if (batch_size > 1) {
    /* HisVal() returns wrong values when the batch size is greater than one,
     * so we calculate the correct values here. */

    scoresum = 0.0;

    if (is_gcu) {
      /* Here, it is assumed that zero batch values are disregarded.
       * This seems to hold for B1_* and DF_* parameters ...
       *
       * Tässä siis oletetaan, että jos batchin arvo on nolla,
       * niin silloin tätä nolla-arvoa ei oteta huomioon ja Nbatch--.
       * Päteekö tämä kaikille gcu parametreille vai ainoastaan vain
       * tietyille (kuten B1_ABS)? */

      Ni = Nbatch;
      j = 0;
      for (i = 0; i < Ni; i++) {
        if (xactive[(i+1)*batch_size - 1] != 0.0) {
          xbatch[j] = xactive[(i+1)*batch_size - 1] - scoresum;
          scoresum += xbatch[j];
          j++;
        }
        else {
          /* Zero batch values are disregarded so we correct Nbatch */
          Nbatch--;
        }
      }
    }
    else {
      /* Not a gcu parameter, zero values are included */

      for (i = 0; i < Nbatch; i++) {
        xbatch[i] = xactive[(i+1)*batch_size - 1] - scoresum;
        scoresum += xbatch[i];
      }
    }
  }
  else {
    for (i = 0; i < Nbatch; i++)
      xbatch[i] = xactive[i];
  }

  /* Nbatch can be zero if, e.g., the B1 iterations didn't converge. The tests
   * are nevertheless performed and they check that Nbatch is positive. */

  /***************************************************************************/


  /***** Perform the statistical tests ***************************************/


  /* Perform the Shapiro-Wilk test */
  
  retsw = ShapiroWilkTest(xbatch, Nbatch, &Z_W_R, &p_W_R);

  /* Perform the D'Agostino-Pearson's K2 test */
  
  retdp = DAgostinoPearsonTest(xbatch, Nbatch, &K2, &p_K2);

  /* Perform the drift-in-mean test (xbatch is divided in half and the
   * means of the halves are compared) */
  
  retdrift2 = DriftInMeanTest(xbatch, Nbatch, Z_drift2, p_drift2, 2);
  
  /* Print the test statistics, p-values and pass values of the tests */
  
  PrintStattests(fp, retsw, Z_W_R, p_W_R, p_W_R >= alpha_sw);
  
  PrintStattests(fp, retdp, K2, p_K2, p_K2 >= alpha_dp);
  
  PrintStattests(fp, retdrift2, Z_drift2[0], p_drift2[0], p_drift2[0] >= alpha_drift);

  fprintf(fp, "\n");

  /***************************************************************************/

  /* Free memory */

  Mem(MEM_FREE, xbatch);
  
}
/*****************************************************************************/


/*****************************************************************************/
static void PrintStattests(FILE *fp, long performed, double Z, double pvalue, long passvalue)
{
  /* Prints the test statistic, the p-value and the pass value of a test 
   * to a file fp. If the test couldn't be performed (performed == 0 or 
   * pvalue == -1), NaN-values are printed instead.
   * 
   * */
  
  if (performed == 0 || pvalue == -1) {
    fprintf(fp, "%16s", "NaN");
    fprintf(fp, "%14s", "NaN");
    fprintf(fp, "%4s", "NaN");
  }
  
  else if (performed == 1) {
    
    /* Print the test statistic */
    
    fprintf(fp, "%16.5E", Z);
    
    /* Print the p-value */
    
    fprintf(fp, "%14.5E", pvalue);

    /* Print the pass value */
    
    fprintf(fp, "%4ld", passvalue);
    
  }
}
/*****************************************************************************/


/*****************************************************************************/
static long ShapiroWilkTest(double *x, long N, double *Z_W_R, double *p_W_R)
{
  /* Shapiro-Wilk normality test with Royston's extension.
   * 
   * The test is only one-tailed (the test passes if p_W_R > alpha, where alpha
   * is the significance level).
   * 
   * Returns 1 if the test was performed, 0 otherwise.
   * */

  double a_1, a_2, a_n, a_nmin1, a_iplus1, axc, S2, W_R, xmean, u, sqrtphi,
      gamma, mu, sigma, mtildess;
  double Nd = (double) N;
  long i, loopstart;
  double *xc, *mtilde;

  /* Polynomial coefficients */
  
  double apoly1[6] = {0.0, 0.221157, -0.147981, -2.07119, 4.434685, -2.706056};
  double apoly2[6] = {0.0, 0.042981, -0.293762, -1.752461, 5.682633, -3.582633};
  double gammapoly1[2] = {-2.273, 0.459};
  double mupoly1[4] = {0.544, -0.39978, 0.025054, -6.714e-4};
  double sigmapoly1[4] = {1.3822, -0.77857, 0.062767, -0.0020322};
  double mupoly2[4] = {-1.5861, -0.31082, -0.083751, 0.0038915};
  double sigmapoly2[3] = {-0.4803, -0.082676, 0.0030302};

  /* The test is valid for N >= 4 */
  
  if(N < 4)
    return 0;

  /* Allocate memory */
  
  xc = (double *)Mem(MEM_ALLOC, N, sizeof(double));
  mtilde = (double *)Mem(MEM_ALLOC, N, sizeof(double));

  /* Copy x to xc */
  
  for (i = 0; i < N; i++)
    xc[i] = x[i];

  /* Sort xc */
  
  qsort(xc, N, sizeof(double), (int(*)(const void*, const void*)) CompareElements);

  /* Return 0 if all the elements of x have the same value */
  
  if (xc[0] == xc[N-1]) {
    Mem(MEM_FREE, xc);
    Mem(MEM_FREE, mtilde);
    return 0;
  }

  /* Calculate the Blom scores and their sum of squares.*/
  
  mtildess = 0.0;
  for (i = 0; i < N; i++) {
    
    if (!NormCDFInv(((double)i + 1.0 - 0.375)/((double)N + 0.25), mtilde + i)) {
      /* NormCDFInv shouldn't fail if everything is done correctly */
      Mem(MEM_FREE, xc);
      Mem(MEM_FREE, mtilde);
      return 0;
    }
    
    mtildess += mtilde[i]*mtilde[i];
  }

  /* Calculate two polynomial coefficients */
  
  apoly1[0] = mtilde[N-1]/sqrt(mtildess);
  apoly2[0] = mtilde[N-2]/sqrt(mtildess);

  /* Calculate the last two and first two elements of vector a */
  
  u = 1.0/sqrt(Nd);
  a_n = Polyval(apoly1, 5, u);
  a_nmin1 = Polyval(apoly2, 5, u);
  a_1 = -a_n;
  a_2 = -a_nmin1;

  /* Calculate the first and last term of the sums axc and S2 */
  
  xmean = SampleMean(xc, N);
  axc = a_1*xc[0] + a_n*xc[N-1];
  S2 = (xc[0] - xmean)*(xc[0] - xmean) + (xc[N-1] - xmean)*(xc[N-1] - xmean);

  /* Calculate sqrt(phi) and set the right start index for the next loop */
  
  if (N <= 5) {
    sqrtphi = sqrt((mtildess - 2.0*mtilde[N-1]*mtilde[N-1])/(1.0 - 2*a_n*a_n));
    loopstart = 1;
  }
  else {
    sqrtphi = sqrt((mtildess - 2.0*mtilde[N-1]*mtilde[N-1]
              - 2.0*mtilde[N-2]*mtilde[N-2])
              /(1.0 - 2.0*a_n*a_n - 2.0*a_nmin1*a_nmin1));
    loopstart = 2;

    /* Calculate the second and the second last term of the sums axc and S2 */
    
    axc += a_2*xc[1] + a_nmin1*xc[N-2];
    S2 += (xc[1] - xmean)*(xc[1] - xmean) + (xc[N-2] - xmean)*(xc[N-2] - xmean);
  }

  /* Calculate the rest of the terms of the sums axc and S2 */
  
  for (i = loopstart; i < N-loopstart; i++) {
    a_iplus1 = mtilde[i]/sqrtphi;
    axc += a_iplus1*xc[i];
    S2 += (xc[i] - xmean)*(xc[i] - xmean);
  }

  /* Calculate the W_R test statistic */
  
  W_R = axc*axc/S2;

  /* Normalize the W_R test statistic using Royston's extension */
  
  if ((N >= 4) && (N <= 11)) {
    gamma = Polyval(gammapoly1, 1, Nd);
    mu = Polyval(mupoly1, 3, Nd);
    sigma = exp(Polyval(sigmapoly1, 3, Nd));
    *Z_W_R = (-log(gamma - log(1.0 - W_R)) - mu)/sigma;
  }
  else {
    mu = Polyval(mupoly2, 3, log(Nd));
    sigma = exp(Polyval(sigmapoly2, 2, log(Nd)));
    *Z_W_R = (log(1.0 - W_R) - mu)/sigma;
  }

  /* Calculate the p-value */
  
  *p_W_R = 1.0 - NormCDF(*Z_W_R);

  /* Free allocated memory */
  Mem(MEM_FREE, xc);
  Mem(MEM_FREE, mtilde);

  return 1;
}
/*****************************************************************************/


/*****************************************************************************/
static long DAgostinoPearsonTest(double *x, long N, double *K2, double *p_K2)
{
  /* D'Agostino-Pearson test for detecting non-normality of the data.
   * The test combines the skewness and kurtosis test statistics given
   * by the functions SkewnessTest() and KurtosisTest().
   * 
   * Returns 1 if the test was performed, 0 otherwise.
   * */

  double Z_b1sqrt, Z_b2, p_tmp, cdf_chi2_2;

  if (SkewnessTest(x, N, &Z_b1sqrt, &p_tmp) 
      && KurtosisTest(x, N, &Z_b2, &p_tmp)) {
    
    /* Test statistic */
    *K2 = Z_b1sqrt*Z_b1sqrt + Z_b2*Z_b2;

    /* K2 has approximately a chi-squared distribution with 2 degrees of 
     * freedom. The exact cumulative distribution function for this case is: */
    cdf_chi2_2 = 1.0 - exp(-*K2/2.0);

    /* p-value of K2 */
    *p_K2 = 1.0 - cdf_chi2_2;

    return 1;
  }
  else {
    return 0;
  }

}
/*****************************************************************************/


/*****************************************************************************/
static long SkewnessTest(double *x, long N, double *Z_b1sqrt, double *p_b1sqrt)
{
  /* Skewness test for detecting non-normality of sample 'x'.
   * 
   * For more information, see:
   * Ralph B. D'agostino, Albert Belanger, Ralph B. D'agostino, "A suggestion 
   * for using powerful and informative tests of normality", 
   * The American Statistician 44, no. 4, 316-321 (1990).
   * 
   * The calculated p-value is two-tailed.
   * 
   * Returns 1 if the test was performed, 0 otherwise.
   * */

  const double Nd = (double)N;
  double m_3, m_2, xmean, b1sqrt, Y, tau, omega2, delta, nu;
  long i;

  /* The test is valid for N >= 9 */
  
  if (N < 9)
    return 0;

  /* Calculate the sample skewness b1sqrt */
  
  xmean = SampleMean(x, N);
  m_3 = 0.0;
  m_2 = 0.0;
  for (i = 0; i < N; i++) {
    m_3 += pow(x[i] - xmean, 3.0);
    m_2 += pow(x[i] - xmean, 2.0);
  }

  m_2 = pow(m_2/Nd, 1.5);

  /* Avoid division by zero (m_2 == 0, when all the elements of x have the 
   * same value) */
  
  if (m_2 == 0.0)
    return 0;
  
  b1sqrt = (m_3/Nd)/m_2;


  /* Parameters for the test statistic */
  
  Y = b1sqrt*sqrt((Nd + 1.0)*(Nd + 3.0)/(6.0*(Nd - 2.0)));
  tau = 3.0*(Nd*Nd + 27.0*Nd - 70.0)*(Nd + 1.0)*(Nd + 3.0) 
        / ((Nd - 2.0)*(Nd + 5.0)*(Nd + 7.0)*(Nd + 9.0));
  omega2 = sqrt(2.0*(tau - 1.0)) - 1;
  
  /* When N is very large (N > 10^15), omega2 -> 1, and the test can't be 
   * performed */
  if (omega2 - 1.0 < ZERO)
    return 0;
  
  delta = 1.0/sqrt(log(sqrt(omega2)));
  nu = sqrt(2.0/(omega2 - 1.0));

  /* Calculate the test statistic and the p-value*/
  
  *Z_b1sqrt = delta*log(Y/nu + sqrt(Y*Y/(nu*nu) + 1.0));
  *p_b1sqrt = 2.0*(1.0 - NormCDF(fabs(*Z_b1sqrt)));

  return 1;
}
/*****************************************************************************/


/*****************************************************************************/
static long KurtosisTest(double *x, long N, double *Z_b2, double *p_b2)
{
  /* Kurtosis test for detecting non-normality of sample 'x'.
   * 
   * For more information, see:
   * Ralph B. D'agostino, Albert Belanger, Ralph B. D'agostino, "A suggestion
   * for using powerful and informative tests of normality", The American
   * Statistician 44, no. 4, 316-321 (1990).
   * 
   * The calculated p-value is two-tailed.
   * 
   * Returns 1 if the test was performed, 0 otherwise.
   * */

  const double Nd = (double)N;
  double m_4, m_2, xmean, b2, Eb2, varb2, xi, beta1sqrt, A;
  long i;

  /* The test is valid for N >= 20 */
  
  if (N < 20)
    return 0;

  /* Calculate the sample kurtosis b2 */
  
  xmean = SampleMean(x, N);
  m_2 = 0.0;
  m_4 = 0.0;
  for (i = 0; i < N; i++) {
    m_4 += pow(x[i] - xmean, 4.0);
    m_2 += pow(x[i] - xmean, 2.0);
  }

  m_2 = pow(m_2/Nd, 2.0);

  /* Avoid division by zero (m_2 == 0, when all the elements of x have the 
   * same value) */
  
  if (m_2 == 0.0)
    return 0;

  b2 = (m_4/Nd)/m_2;


  /* Parameters for the test statistic */
  
  Eb2 = 3.0*(Nd - 1.0)/(Nd + 1.0);
  varb2 = 24.0*Nd*(Nd - 2.0)*(Nd - 3.0)
          /((Nd + 1.0)*(Nd + 1.0)*(Nd + 3.0)*(Nd + 5.0));
  xi = (b2 - Eb2)/sqrt(varb2);
  beta1sqrt = 6.0*(Nd*Nd - 5.0*Nd + 2.0)/((Nd + 7.0)*(Nd + 9.0))
              *sqrt(6.0*(Nd + 3.0)*(Nd + 5.0) / (Nd*(Nd - 2.0)*(Nd - 3.0)));
  A = 6.0 + 8.0/beta1sqrt*(2.0/beta1sqrt 
      + sqrt(1.0 + 4.0/(beta1sqrt*beta1sqrt)));

  /* Calculate the test statistic and the two-tailed p-value */
  
  if (1.0 + xi*sqrt(2.0/(A - 4.0)) <= ZERO) {
    /* In difficult cases (e.g. x has two peaks far from each other),
     * 1 + xi*sqrt(2/(A - 4)) <= 0 and pow() gives NaN. Therefore we set: */
    *Z_b2 = -INFTY;
    *p_b2 = 0.0;
  }
  else {
    *Z_b2 = sqrt(9.0*A/2.0)*((1.0 - 2.0/(9.0*A)) 
            - pow((1.0 - 2.0/A)/(1.0 + xi*sqrt(2.0/(A - 4.0))), 1.0/3.0));
    *p_b2 = 2.0*(1.0 - NormCDF(fabs(*Z_b2)));
  }

  return 1;
}
/*****************************************************************************/


/*****************************************************************************/
static long DriftInMeanTest(double *x, long N, double *Z_drift, 
                            double *p_drift, long m)
{
  /* Performs the drift-in-mean test (unequal variance t-test).
   * The data vector x is divided into m parts, and the test is
   * performed for the adjacent parts. Each part must have at least two elements.
   * The calculated p_drift is the two-tailed p-value.
   * 
   * Returns 1 if the test was performed, 0 otherwise.
   * If the test fails for some of the compared parts, the corresponding p_drift 
   * is set to -1 or -2. -1 means that all the elements in either of the 
   * compared parts were the same. -2 means that tCDF()- or tCDFNormalApprox()-
   * function failed, which shouldn't happen.
   * 
   * */
  
  double *xi, *mean, *var;
  long *Ni;
  double nu, u, tmp, tcdf, Nid, Nip1d;
  long Nm, modNm, rettcdf, i, j;  /* Nm = number of elements in parts */

  /* The test can't be performed in these cases: */
  
  if ((N < 4) || (m < 2) || (2*m > N))
    return 0;


  mean = (double *)Mem(MEM_ALLOC, m, sizeof(double)); /* Means of the parts */
  var = (double *)Mem(MEM_ALLOC, m, sizeof(double));  /* Unbiased sample variances of the parts */
  Ni = (long *)Mem(MEM_ALLOC, m, sizeof(long));  /* Number of elements in the parts */

  Nm = N/m;
  modNm = N - Nm*m;


  /* Loop over parts*/
  
  xi = x;
  for (i = 0; i < m; i++) {

    /* The number of elements in a part is either Nm or Nm+1. 
     * (Ni[i] > 1, because it was checked that N >= 4 && 2*m <= N) */
    
    if (i < modNm)
      Ni[i] = 1 + Nm;
    else
      Ni[i] = Nm;
    
    /* Calculate the mean */
    
    mean[i] = SampleMean(xi, Ni[i]);
    
    /* Calculate the variance */
    
    var[i] = 0.0;
    for (j = 0; j < Ni[i]; j++)
      var[i] += (xi[j] - mean[i])*(xi[j] - mean[i]);
    var[i] /= ((double)Ni[i] - 1.0);

    xi += Ni[i];
  }

  /* Calculate the test statistics and the p-values */
  
  for (i = 0; i < m-1; i++) {
    
    if ((var[i] == 0.0) || (var[i+1] == 0.0)) {
      /* var[i] == 0 when all the elements in xi are the same */
      p_drift[i] = -1;
    }
    
    else {
      Nid = (double)Ni[i];
      Nip1d = (double)Ni[i+1];
      
      Z_drift[i] = (mean[i+1] - mean[i])/sqrt(var[i]/Nid + var[i+1]/Nip1d);
      
      /* Z_drift is t-distributed with nu degrees of freedom */
      
      u = var[i+1]/var[i];
      tmp = 1.0/Nid + u/Nip1d;
      nu = tmp*tmp / (1.0/(Nid*Nid*(Nid - 1.0)) + 
                      u*u/(Nip1d*Nip1d*(Nip1d - 1.0)));
      
      if (nu < 400) {
        /* Use the 'accurate' cdf of Student's t-distribution. */
        rettcdf = tCDF(fabs(Z_drift[i]), nu, &tcdf);
      }
      else {
        /* Use the Gaver's and Kafadar's normal approximation. For nu = 400, 
         * the relative error for calculating the p-value 0.01 is
         * around 3.5 × 10−6 */
        rettcdf = tCDFNormalApprox(fabs(Z_drift[i]), nu, &tcdf);
      }
      
      
      /* Calculate the two-tailed p-value */
      
      if (rettcdf == 0) {
        /* Something went wrong in tCDF() or tCDFNormalApprox(). 
         * Normally this should not happen. */
        p_drift[i] = -2;
      }
      else {
        p_drift[i] = 2.0*(1.0 - tcdf);
      }
      
    }
  }

  Mem(MEM_FREE, mean);
  Mem(MEM_FREE, var);
  Mem(MEM_FREE, Ni);
  
  return 1;

}
/*****************************************************************************/


/*****************************************************************************/
static long tCDF(double x, double nu, double *cdf)
{
  /* Calculates the cumulative distribution function (cdf) of Student's 
   * t-distribution. The cdf is calculated using the regularized incomplete 
   * beta function I_x(a,b).
   * 
   * With this implementation, the relative error (RE) of the calculated
   * cdf is 
   * RE < 1.1*10^-7, 
   * when 0.5 < nu < 400 and -50 < x < 50.
   * 
   * Returns 1 if successful, 0 otherwise.
   * 
   * THIS FUNCTION HAS BEEN TESTED ONLY AT THE INTERVALS 
   * 0.5 < nu < 400 and -50 < x < 50.
   * */
  
  double a, b, nux, tol;
  long maxit;
  
  /* Check that nu is positive */
  
  if (nu <= 0.0)
    return 0;
  
  /* Parameteres for Betareg */
  
  a = nu/2.0;
  b = 0.5;
  nux = nu/(x*x + nu);
  tol = 1.0e-10/a;        /* Tolerance for betareg */
  maxit = 1000;           /* Maximum number of terms summed in betareg */


  if (nux < a/(a + b + 2.0)) {
    /* Calculate I_x(a,b). (Matlab uses nux < (a+1)/(a+b+2). However,
     * nux < a/(a+b+2) produces more accurate result in this case.) */
      
    if (Betareg(nux, a, b, tol, maxit, cdf) == 0) {
      /* Betareg failed. Doesn't happen when 0.5 < nu < 400 and 
        * -50 < x < 50. */
      return 0;
    }
  }
  else {
    /* Calculate I_x(a,b) = 1 - I_{1-x}(b,a) */
    
    if (Betareg(1.0 - nux, b, a, tol, maxit, cdf) == 0) {
      /* Betareg failed. Doesn't happen when 0.5 < nu < 400 and 
        * -50 < x < 50. */
      return 0;
    }
    *cdf = 1.0 - *cdf;
  }

  *cdf /= 2.0;

  if (x > 0.0) {
      *cdf = 1.0 - *cdf;
  }
  else if (x == 0.0) {
      /* The accurate result */
      *cdf = 0.5;
  }

  return 1;
}
/*****************************************************************************/


/*****************************************************************************/
static long tCDFNormalApprox(double x, double nu, double *cdf)
{
  /* Calculates the cumulative distribution function (cdf) of Student's 
   * t-distribution using the Gaver-Kafadar normal approximation 
   * (1984, Amer. Statist. 38, 308–311).
   * 
   * Returns 1 if successful, 0 otherwise.
   * */
  
  double g, z;
  
  /* Check that nu is positive */
  
  if (nu <= 0.0)
    return 0;

  /* Gaver-Kafadar normal approximation */

  g = (nu - 1.5)/((nu - 1.0)*(nu - 1.0));
  z = sqrt(log(1.0 + x*x/nu)/g);

  if (x < 0.0)
    *cdf = NormCDF(-z);
  else
    *cdf = NormCDF(z);

  return 1;
}
/*****************************************************************************/


/*****************************************************************************/
static long Betareg(double x, double a, double b, double tol, long maxit, 
                    double *betareg)
{
  /* Calculates the regularized incomplete beta function I_x(a,b) which is the
   * ratio of the incomplete and complete beta functions. The incomplete beta
   * function is calculated by truncating its series after maxit terms or when 
   * the tolerance tol is reached.
   * 
   * Returns 1 if successful, 0 otherwise.
   * 
   * THIS FUNCTION HAS BEEN TESTED ONLY WITH THE FUNCTION tCDF().
   * */
  
  double tmp, sum, betacomp, betaincomp;
  long j;
  
  /* Check the input values */
  
  if (x < 0.0 || x > 1.0 || a <= 0.0 || b <= 0.0)
    return 0;

  tmp = 1.0;
  sum = 0.0;
  j = 0;
  
  /* Calculate the truncated series */
  
  while ((fabs(tmp) > tol) && (j < maxit)) {
    j = j + 1;
    tmp = tmp*(1.0 - b/(double)j)*x;
    sum = sum + tmp/(a + (double)j);
  }

  /* Complete beta function */
  
  betacomp = exp(lgamma(a) + lgamma(b) - lgamma(a + b));

  /* Incomplete beta function */
  
  betaincomp = pow(x, a)/a*(1.0 + a*sum);
  
  *betareg = betaincomp/betacomp;
  
  /* If the convergence of the series is slow, betareg can be above one or 
   * below zero. It has been tested that this doesn't happen with the 
   * implemented version of tCDF() when nu < 400. */
  
  if (*betareg > 1.0 || *betareg < 0.0)
    return 0;

  return 1;
}
/*****************************************************************************/


/*****************************************************************************/
static double SampleMean(double *x, long N)
{
  /* Returns the sample mean of x.
   * */
  
  double xmean = 0.0;
  long i;

  for (i = 0; i < N; i++)
    xmean += x[i];
  return xmean/(double)N;
}
/*****************************************************************************/


/*****************************************************************************/
static int CompareElements(const double *a, const double *b)
{
  /* Compare elements *a and *b. 
   * 
   * Returns 1 if a > b, 0 if a == b, or -1 if a < b.
   * */
  
  double diff = *a-*b;
  if (diff < 0)
    return -1;
  else if (diff > 0)
    return 1;
  else
    return 0;
}
/*****************************************************************************/


/*****************************************************************************/
static double Polyval(const double *coeff, long degree, double x)
{
  /* Evaluates a polynomial of degree 'degree' with coefficients 'coeff'
   * at 'x'.
   * */
  
  double p, f;
  long i;

  f = coeff[0];
  if (degree > 0) {
    p = x*coeff[degree];
    for (i = degree-1; i > 0; i--)
      p = (p + coeff[i])*x;
    f += p;
  }
  return f;
}
/*****************************************************************************/


/*****************************************************************************/
static double NormCDF(double x)
{
  /* Returns the value of the cumulative distribution function of the
   * standard normal distribution.
   * */
  
  return 0.5*(1.0 + erf(x/sqrt(2.0)));
}
/*****************************************************************************/


/*****************************************************************************/
static long NormCDFInv(double p, double *Z)
{
  /* Returns the inverse of the cumulative distribution function
   * of the standard normal distribution. The used algorithm is:
   * 
   * ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3, 477- 484
   * 
   * which has the description:
   * "Produces the normal deviate Z corresponding to a given lower
   * tail area of P; Z is accurate to about 1 part in 10**16."
   * 
   * Original Fortran code: http://jblevins.org/mirror/amiller/as241.f90
   * 
   * Returns 1 if successful, 0 otherwise.
   * 
   * */
  
  double q, r;
  double split1 = 0.425, split2 = 5.0, const1 = 0.180625, const2 = 1.6;

  /* Coefficients for p close to 0.5 */
  
  double a0 = 3.3871328727963666080e0,
         a1 = 1.3314166789178437745e+2,
         a2 = 1.9715909503065514427e+3,
         a3 = 1.3731693765509461125e+4,
         a4 = 4.5921953931549871457e+4,
         a5 = 6.7265770927008700853e+4,
         a6 = 3.3430575583588128105e+4,
         a7 = 2.5090809287301226727e+3,
         b1 = 4.2313330701600911252e+1,
         b2 = 6.8718700749205790830e+2,
         b3 = 5.3941960214247511077e+3,
         b4 = 2.1213794301586595867e+4,
         b5 = 3.9307895800092710610e+4,
         b6 = 2.8729085735721942674e+4,
         b7 = 5.2264952788528545610e+3;

  /* Coefficients for p not close to 0, 0.5 or 1 */
  
  double c0 = 1.42343711074968357734e0,
         c1 = 4.63033784615654529590e0,
         c2 = 5.76949722146069140550e0,
         c3 = 3.64784832476320460504e0,
         c4 = 1.27045825245236838258e0,
         c5 = 2.41780725177450611770e-1,
         c6 = 2.27238449892691845833e-2,
         c7 = 7.74545014278341407640e-4,
         d1 = 2.05319162663775882187e0,
         d2 = 1.67638483018380384940e0,
         d3 = 6.89767334985100004550e-1,
         d4 = 1.48103976427480074590e-1,
         d5 = 1.51986665636164571966e-2,
         d6 = 5.47593808499534494600e-4,
         d7 = 1.05075007164441684324e-9;

  /* Coefficients for p near 0 or 1. */
  
  double e0 = 6.65790464350110377720e0,
         e1 = 5.46378491116411436990e0,
         e2 = 1.78482653991729133580e0,
         e3 = 2.96560571828504891230e-1,
         e4 = 2.65321895265761230930e-2,
         e5 = 1.24266094738807843860e-3,
         e6 = 2.71155556874348757815e-5,
         e7 = 2.01033439929228813265e-7,
         f1 = 5.99832206555887937690e-1,
         f2 = 1.36929880922735805310e-1,
         f3 = 1.48753612908506148525e-2,
         f4 = 7.86869131145613259100e-4,
         f5 = 1.84631831751005468180e-5,
         f6 = 1.42151175831644588870e-7,
         f7 = 2.04426310338993978564e-15;


  /* Check the boundaries of p */
  
  if ((p < 0.0) || (p > 1.0)) {
    return 0;
  }
  else if (p < ZERO) {
    *Z = -INFTY;
    return 1;
  }
  else if (p == 1.0) {
    *Z = INFTY;
    return 1;
  }

  q = p - 0.5;

  if (fabs(q) < split1) {
    /* 0.075 < p < 0.925 */
    
    r = const1 - q*q;
    *Z = q*(((((((a7*r + a6)*r + a5)*r + a4)*r + a3)*r + a2)*r + a1)*r + a0)
         /(((((((b7*r + b6)*r + b5)*r + b4)*r + b3)*r + b2)*r + b1)*r + 1.0);
  }
  else {
    /* min(p,1-p) <= 0.075 */
    
    if (q < 0.0)
      r = p;
    else
      r = 1.0 - p;
    
    r = sqrt(-log(r)); /* min(p,1-p) = exp(-r^2) */

    if (r <= split2) {
      /* min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
      r -= const2;
      *Z = (((((((c7*r + c6)*r + c5)*r + c4)*r + c3)*r + c2)*r + c1)*r + c0)
           /(((((((d7*r + d6)*r + d5)*r + d4)*r + d3)*r + d2)*r + d1)*r + 1.0);
    }
    else {
      /* min(p,1-p) < exp(-25) */
      r -= split2;
      *Z = (((((((e7*r + e6)*r + e5)*r + e4)*r + e3)*r + e2)*r + e1)*r + e0)
           /(((((((f7*r + f6)*r + f5)*r + f4)*r + f3)*r + f2)*r + f1)*r + 1.0);
    }
    
    if (q < 0.0)
      *Z = -*Z;
  }
  
  return 1;
}
/*****************************************************************************/
