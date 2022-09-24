/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : mcdf.c                                         */
/*                                                                           */
/* Created:       2013/03/08 (JLe)                                           */
/* Last modified: 2013/03/08 (JLe)                                           */
/* Version:       2.1.13                                                     */
/*                                                                           */
/* Description: Laskee epäjatkuvuustekijöiden laskentaan tarvittavan         */
/*              homogeenisen vuon moniryhmäisellä Monte Carlolla             */
/*                                                                           */
/* Comments: - Onko tolla reunavirralla paikkariippuvuus, vai                */
/*             homogenisoidaanko se pois?                                    */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MCDF:"

#define POLTTOAINE 1
#define REFLEKTORI 2

void Sironta(double *, double *, double *, long *, double *, long, double, 
	     double, const double *, const double *, long);

void Lahde(double *, double *, double *, double *, double *, double *, long *, 
	   double *, double *, double *, long, long, long, const double *, 
	   long, long);

void Lahde2(double *, double *, double *, double *, double *, double *,
	    long *, double *, double **, long);

double ArvoP1Mu(double, long);

long RegIdx(double x, double y, double z, double u, double v, double w,
	    long *cell, long id);

void Sironta2(double *u, double *v, double *w, double *wgt, long *n0, long nmg,
	      double mubar, double *scattp, double *scattw, long id);

void Fissio(double x, double y, double z, double wgt, double nubar, 
	    double *chi, long id);

/*****************************************************************************/

void MCDF()
{
#ifdef mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm

  long gcu, adf, ptr, nmg, surf, ns, type, mode, n, m, ng, id, n1, n2, m1, m2;
  long pop, nbatch, i, uni, ngu, ncyc, skip;
  const double *abs, *fiss, *nsf, *smtx, *pmtx, *mubar, *Iout, *het;
  const double *params;
  double *hom, *scatt, *tot, *alb, *src, *mu0, sum, x, y, z, u, v, w, wgt, l;
  double d, f, u0, v0, w0, mu, norm, *Iin, E, t;
  static double dis[70][100];
  FILE *fp;
  char tmpstr[MAX_STR];
  double *soos[2];

  return;

#ifdef mmmmmmmmmmmmmmmmmmmmmmm

  double **utot, **uabs, **ufiss, **uscatt, **unubar, **uchi, **umubar;
  double ***uscattp, ***uscattw, *majorant;
  long cell, nn;

  /* Check active cycle */
  
  if (RDB[DATA_CYCLE_IDX] < RDB[DATA_CRIT_SKIP])
    return;

  /* Check if group constants are generated */

  if((long)RDB[DATA_OPTI_GC_CALC] == NO)
    return;

  /* Check micro-group mode */

  if (((long)RDB[DATA_SIMULATION_COMPLETED] == NO) && 
      ((long)RDB[DATA_MICRO_CALC_MODE] == MICRO_CALC_MODE_END))
    return;
  else if (((long)RDB[DATA_SIMULATION_COMPLETED] == YES) && 
	   ((long)RDB[DATA_MICRO_CALC_MODE] == MICRO_CALC_MODE_CYCLE))
    return;

  /***************************************************************************/

  /***** Calculate cross sections ********************************************/

  /* Reset universe indexes */

  uni = (long)RDB[DATA_PTR_U0];
  while (uni > VALID_PTR)
    {
      /* Reset index */

      WDB[uni + UNIVERSE_GCU_IDX] = -1.0;
      
      /* Next universe */

      uni = NextItem(uni);
    }

  /* Get number of gc universes */

  gcu = (long)RDB[DATA_PTR_GCU0];
  ngu = ListSize(gcu);

  /* Allocate memory for lists */

  utot = (double **)Mem(MEM_ALLOC, ngu, sizeof(double *));  
  uabs = (double **)Mem(MEM_ALLOC, ngu, sizeof(double *));  
  ufiss = (double **)Mem(MEM_ALLOC, ngu, sizeof(double *));  
  uscatt = (double **)Mem(MEM_ALLOC, ngu, sizeof(double *));  
  unubar = (double **)Mem(MEM_ALLOC, ngu, sizeof(double *));  
  uchi = (double **)Mem(MEM_ALLOC, ngu, sizeof(double *));  
  umubar = (double **)Mem(MEM_ALLOC, ngu, sizeof(double *));  
  uscattp = (double ***)Mem(MEM_ALLOC, ngu, sizeof(double **));  
  uscattw = (double ***)Mem(MEM_ALLOC, ngu, sizeof(double **));  

  /* Get number of energy groups */
	      
  nmg = (long)RDB[DATA_ERG_FG_NG];

  /* Allocate memory for majorant */

  majorant = (double *)Mem(MEM_ALLOC, nmg, sizeof(double));
  memset(majorant, 0.0, nmg*sizeof(double));

  /* Loop over universes */
  
  i = 0;
  gcu = (long)RDB[DATA_PTR_GCU0];
  while (gcu > VALID_PTR)
    {
      /* Allocate memory for cross sections */

      utot[i] = (double *)Mem(MEM_ALLOC, nmg, sizeof(double));  
      uabs[i] = (double *)Mem(MEM_ALLOC, nmg, sizeof(double));  
      ufiss[i] = (double *)Mem(MEM_ALLOC, nmg, sizeof(double));  
      uscatt[i] = (double *)Mem(MEM_ALLOC, nmg, sizeof(double));  
      unubar[i] = (double *)Mem(MEM_ALLOC, nmg, sizeof(double));  
      uchi[i] = (double *)Mem(MEM_ALLOC, nmg, sizeof(double));  
      umubar[i] = (double *)Mem(MEM_ALLOC, nmg, sizeof(double));  

      /* Allocate memory for scattering matrixes */

      uscattp[i] = (double **)Mem(MEM_ALLOC, nmg, sizeof(double *));  
      for (n = 0; n < nmg; n++)
	uscattp[i][n] = (double *)Mem(MEM_ALLOC, nmg, sizeof(double));  

      uscattw[i] = (double **)Mem(MEM_ALLOC, nmg, sizeof(double *));  
      for (n = 0; n < nmg; n++)
	uscattw[i][n] = (double *)Mem(MEM_ALLOC, nmg, sizeof(double));  

      /* Put data */

      for (n = 0; n < nmg; n++)
	{
	  /* Directly available values */

	  ptr = (long)RDB[gcu + GCU_RES_FG_ABSXS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  uabs[i][n] = Mean(ptr, n + 1);

	  ptr = (long)RDB[gcu + GCU_RES_FG_FISSXS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  ufiss[i][n] = Mean(ptr, n + 1);

	  ptr = (long)RDB[gcu + GCU_RES_FG_SCATTXS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  uscatt[i][n] = Mean(ptr, n + 1);

	  ptr = (long)RDB[gcu + GCU_RES_FG_NUBAR];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  unubar[i][n] = Mean(ptr, n + 1);

	  ptr = (long)RDB[gcu + GCU_RES_FG_CHI];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  uchi[i][n] = Mean(ptr, n);
	  
	  ptr = (long)RDB[gcu + GCU_RES_FG_P1_MUBAR];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  umubar[i][n] = Mean(ptr, n + 1);
	  
	  /* Calculate total */

	  utot[i][n] = uabs[i][n] + ufiss[i][n] + uscatt[i][n];

	  /* Compare to majorant */

	  if (utot[i][n] > majorant[n])
	    majorant[n] = utot[i][n];

	  /* Group transfer probability matrix */

	  ptr = (long)RDB[gcu + GCU_RES_FG_GTRANSP];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  
	  for (m = 0; m < nmg; m++)
	    uscattp[i][n][m] = Mean(ptr, m, n);

	  /* Group transfer weight matrix */

	  ptr = (long)RDB[gcu + GCU_RES_FG_GPRODP];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  
	  for (m = 0; m < nmg; m++)
	    {
	      if (uscattp[i][n][m] > 0.0)
		uscattw[i][n][m] = Mean(ptr, m, n)/uscattp[i][n][m];
	      else
		uscattw[i][n][m] = 0.0;
	    }
	}
   
      /* Pointer to universe */

      uni = (long)RDB[gcu + GCU_PTR_UNIV];
      CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

      printf("\nUniverse %s:\n", GetText(uni + UNIVERSE_PTR_NAME));

      for (n = 0; n < nmg; n++)
	{
	  printf("%E %E %E %E\n", utot[i][n], uabs[i][n], ufiss[i][n], 
		 uscatt[i][n]);
	  printf("%E %E %E\n", uchi[i][n], unubar[i][n], umubar[i][n]); 

	  for (m = 0; m < nmg; m++)
	    printf("%E ", uscattp[i][n][m]);
	  printf("\n");

	  for (m = 0; m < nmg; m++)
	    printf("%E ", uscattw[i][n][m]);
	  printf("\n");
	}

      /* Put index */

      WDB[uni + UNIVERSE_GCU_IDX] = (double)(i++);
      
      /* Next universe */

      gcu = NextItem(gcu);
    }

  /***************************************************************************/

  /***** Generate initial source *********************************************/

  /* Reset OpenMP thread id */

  id = 0;

  ptr = (long)RDB[DATA_PART_PTR_SOURCE];
  printf("source = %ld\n", ListSize(ptr) - 1);
  ptr = (long)RDB[OMPPtr(DATA_PART_PTR_NSTACK, id)];
  printf("stack = %ld\n", ListSize(ptr) - 1);
  ptr = (long)RDB[OMPPtr(DATA_PART_PTR_BANK, id)];
  printf("bank = %ld\n", ListSize(ptr) - 1);
  ptr = (long)RDB[OMPPtr(DATA_PART_PTR_QUE, id)];
  printf("que = %ld\n", ListSize(ptr) - 1);

  /* Flush bank */
  /*
  FlushBank();
  */
  /* Get population size */

  pop = 1000;

  /* Loop over population */

  for (n = 0; n < pop; n++)
    {
      /* Re-sampling loop */

      for (m = 0; m < 1000000; m++)
	{
	  /* Sample coordinates */
      
	  x = RandF(0)*(RDB[DATA_GEOM_MAXX] - RDB[DATA_GEOM_MINX]) + 
	    RDB[DATA_GEOM_MINX];
	  
	  y = RandF(0)*(RDB[DATA_GEOM_MAXY] - RDB[DATA_GEOM_MINY]) + 
	    RDB[DATA_GEOM_MINY];
	  
	  if ((long)RDB[DATA_GEOM_DIM] == 3)
	    {
	      z = RandF(0)*(RDB[DATA_GEOM_MAXZ] - RDB[DATA_GEOM_MINZ]) + 
		RDB[DATA_GEOM_MINZ];
	    }
	  else
	    z = 0.0;

	  /* Find region */

	  i = RegIdx(x, y, z, u, v, w, &cell, id);
	  
	  /* Check fissile */
	  
	  if (i > 0)
	    if (ufiss[i][0] > 0.0)
	      break;
	}

      /* Check number of resamples */

      if (m == 1000000)
	Die(FUNCTION_NAME, "Unable to sample initial source");
      
      /* Sample direction cosines */
      
      IsotropicDirection(&u, &v, &w, id);

      /* Get particle from stack */

      nn = FromStack(PARTICLE_TYPE_NEUTRON, id);
  
      /* Put values */
      
      WDB[nn + PARTICLE_X] = x;
      WDB[nn + PARTICLE_Y] = y;
      WDB[nn + PARTICLE_Z] = z;
      
      WDB[nn + PARTICLE_U] = u;
      WDB[nn + PARTICLE_V] = v;
      WDB[nn + PARTICLE_W] = w;
      
      WDB[nn + PARTICLE_E] = 0.0;
      WDB[nn + PARTICLE_WGT] = 1.0;
      WDB[nn + PARTICLE_FMTX_IDX] = -1.0;

      /* Put particle to bank */
      
      ToBank(nn, id);
    }

  /* Put population size */
  
  WDB[DATA_SIMUL_BATCH_SIZE] = (double)pop;

  /***************************************************************************/

  /***** Transport cycle *****************************************************/

  /* Number of active and inactive cycles */

  ncyc = 100;
  skip = 10;

  /* Reset initial guess for k-eff */

  WDB[DATA_CYCLE_KEFF] = 1.0;

  /* Loop over cycles */

  for (n = 0; n < ncyc + skip; n++)
    {
      /* Put cycle index */

      WDB[DATA_CYCLE_IDX] = (double)n;

      /* Clear buffer */
      
      ClearBuf();
      
      /* Normalize source */
      
      NormalizeCritSrc();

      printf("%ld %1.5f\n", n + 1, RDB[DATA_CYCLE_KEFF]);

      /* Loop over source */

      while (FromSrc(id) > VALID_PTR)
	{
	  /* Get particle from que */

	  nn = FromQue(id);

	  /* Get variables */

	  x = RDB[nn + PARTICLE_X];
	  y = RDB[nn + PARTICLE_Y];
	  z = RDB[nn + PARTICLE_Z];
	  u = RDB[nn + PARTICLE_U];
	  v = RDB[nn + PARTICLE_V];
	  w = RDB[nn + PARTICLE_W];
	  ng = (long)RDB[nn + PARTICLE_E];
	  wgt = RDB[nn + PARTICLE_WGT];

	  /* Transport loop */

	  while (1 != 2)
	    {
	      /* Sample distance */

	      d = -log(RandF(id))/majorant[ng];

	      /* Move neutron */

	      x = x + u*d;
	      y = y + v*d;
	      z = z + w*d;

	      /* Find region */
	  
	      i = RegIdx(x, y, z, u, v, w, &cell, id);

	    

	      /* Check leak */

	      if (i < 0)
		{
		  /* Apply boundary conditions */
	      
		  if (BoundaryConditions(&cell, &x, &y, &z, &u, &v, &w, 
					 1.0, &wgt, id) < 0)
		    {
		      printf("%E %E %E\n", x, y, z);

		      /* Point is outside, break loop */

		      break;
		    }
		  
		  /* Cycle loop */

		  continue;
		}

	      /* Sample between real and virtual collision */
	    
	      if (RandF(id) > utot[i][ng]/majorant[ng])
		{
		  /* Virtual collision, cycle loop */

		  continue;
		}

	      /* Sample fraction of total cross section */
			  
	      f = RandF(id)*utot[i][ng];

	      /* Find reaction mode */
	      
	      if ((f = f - uscatt[i][ng]) < 0.0)
		{
		  /* Scattering */

		  Sironta2(&u, &v, &w, &wgt, &ng, nmg, umubar[i][ng], 
			   uscattp[i][ng], uscattw[i][ng], id);
	
		  /* Cycle loop */

		  continue;
		}
	      else if ((f = f - uabs[i][ng]) < 0.0)
		{
		  /* Absorption, break loop */

		  break;
		}
	      else if ((f = f - ufiss[i][ng]) < 0.0)
		{
		  /* Fission */

		  Fissio(x, y, z, wgt, unubar[i][ng], uchi[i], id);

		  /* Break loop */
		  
		  break;		  
		}
	      else
		Die(FUNCTION_NAME, "Unable to sample reaction");
	    }
	  
	  /* History terminated, put particle back in stack */

	  ToStack(nn, id);
	}
    }

  /***************************************************************************/

#endif


  fprintf(out, "Lasketaan epäjatkuvuustekijät moniryhmä Monte Carlolla...\n");

  sprintf(tmpstr, "%s_fg.m", GetText(DATA_PTR_INPUT_FNAME));
  fp = fopen(tmpstr, "w");

  /* Loop over universes */

  gcu = (long)RDB[DATA_PTR_GCU0];
  while (gcu > VALID_PTR)
    {
      /* Haetaan pointteri epäjatkuvuustekijöihin (ja vaihdetaan samalla */
      /* salakieleen kiinalaisten kiusaksi...) */

      if ((adf = (long)RDB[gcu + GCU_PTR_ADF]) < VALID_PTR)
	{
	  /* Next universe */

	  gcu = NextItem(gcu);

	  /* Cycle loop */

	  continue;
	}

      /***********************************************************************/

      /***** Haetaan laskussa tarvittavat moniryhmävakiot ********************/

      /* Pointteri mikroryhmäjakoon */
  
      ptr = (long)RDB[DATA_MICRO_PTR_EGRID];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      
      /* Mikroryhmien lukumäärä */
      
      nmg = (long)RDB[ptr + ENERGY_GRID_NE] - 1;

      /* Rajapintojen lukumäärä */

      ns = (long)RDB[adf + ADF_NSURF];

      /* Pointteri rajapintaan */

      surf = (long)RDB[adf + ADF_PTR_SURF];
      CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

      /* Pinnan tyyppi */

      type = (long)RDB[surf + SURFACE_TYPE];

      /* Parametrivektori */

      ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      params = &RDB[ptr];

      /* Pointterit vaikutusaloihin */

      ptr = (long)RDB[gcu + GCU_MICRO_ABS];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      abs = &RES2[ptr];
      
      ptr = (long)RDB[gcu + GCU_MICRO_FISS];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      fiss = &RES2[ptr];
      
      ptr = (long)RDB[gcu + GCU_MICRO_NSF];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      nsf = &RES2[ptr];
      
      ptr = (long)RDB[gcu + GCU_MICRO_SCATT_MTX];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      smtx = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_SCATT_PROD_MTX];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      pmtx = &RES2[ptr];
      
      ptr = (long)RDB[gcu + GCU_MICRO_MUBAR];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      mubar = &RES2[ptr];
      
      ptr = (long)RDB[gcu + GCU_MICRO_ADF_IN_CURR];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      Iin = &RES2[ptr];
      
      ptr = (long)RDB[gcu + GCU_MICRO_ADF_OUT_CURR];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      Iout = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_ADF_SURF_FLUX];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      het = &RES2[ptr];

      /* Varataan muistia sirontavaikutusalalle albedoille ja homogeeniselle */
      /* vuoratkaisulle (käytetään myöhemmin esivarattua muistia) */

      scatt = (double *)Mem(MEM_ALLOC, nmg, sizeof(double));
      tot = (double *)Mem(MEM_ALLOC, nmg, sizeof(double));
      alb = (double *)Mem(MEM_ALLOC, ns*nmg, sizeof(double));
      hom = (double *)Mem(MEM_ALLOC, ns*nmg, sizeof(double));
      mu0 = (double *)Mem(MEM_ALLOC, ns*nmg, sizeof(double));
      
      /* Tätä ei pidä varata kriittisyysmoodissa */

      src = (double *)Mem(MEM_ALLOC, ns*nmg, sizeof(double));

      /* Asetetaan laskentamoodi */

      mode = REFLEKTORI;

      /* Luuppi energiaryhmien yli */

      for (n = 0; n < nmg; n++)
	{
	  /* Lasketaan sirontavaikutusala */

	  scatt[n] = 0.0;
	  for (m = 0; m < nmg; m++)
	    scatt[n] = scatt[n] + smtx[n*nmg + m];
	  
	  /* Lasketaan kokonaisvaikutusala */

	  tot[n] = abs[n] + fiss[n] + scatt[n];

	  /* Vaihdetaan moodia jos fissiili */

	  if (fiss[n] > 0.0)
	    mode = POLTTOAINE;
	}
      
      /* Lasketaan neutronien keskimääräinen suuntakosini reunalla */

      for (n = 0; n < nmg*ns; n++)
	{
	  if (het[n] > 0.0)
	    mu0[n] = (Iout[n] - Iin[n])/het[n];
	  else
	    mu0[n] = 0.0;
	}

      /* Lasketaan albedot */
	  
      for (n = 0; n < nmg*ns; n++)
	{
	  if (Iout[n] > 0.0)
	    alb[n] = Iin[n]/Iout[n];
	  else
	    alb[n] = 1.0;
	}

      printf("tot, fiss, abs, scatt:\n");
      
      for (n = 0; n < nmg; n++)
	printf("%3ld %1.2E  %1.2E  %1.2E  %1.2E %1.5f %5.2f\n", 
	       n, tot[n], fiss[n], abs[n], scatt[n], 
	       scatt[n]/(tot[n] - abs[n]), mubar[n]);

      printf("in, out, het, alb:\n");
      for (n = 0; n < nmg*ns; n++)
	printf("%3ld %1.2E %1.2E %1.2E %1.2E\n", n, Iin[n], Iout[n], het[n], alb[n]);

      /***********************************************************************/
      
      /***** Ajetaan moniryhmä Monte Carlo -simulaatio ***********************/

      /* Tarkistetaan moodi */

      if (mode == POLTTOAINE)
	{
	  /***** Kriittisyyslähdesimulaatio **********************************/

	  /* Ei toimi vielä */

	  Die(FUNCTION_NAME, "Kriittisyyslähdesimulaatio ei toimi vielä");

	  /*******************************************************************/
	}
      else
	{
	  /***** Ulkoisen lähteen simulaatio *********************************/

	  for (n = 0; n < nmg; n++)
	    for (i = 0; i < 100; i++)
	      dis[n][i] = 0.0;

	  if (1 != 2)
	    {
	      /* Lasketaan kokonaisvirta sisään */

	      sum = 0;
	      for (n = 0; n < nmg*ns; n++)
		sum = sum + Iin[n];
	      
	      /* Tarkistetaan */
	      
	      if (sum < ZERO)
		Die(FUNCTION_NAME, "Jotain häikkää");
	      
	      /* Lasketaan lähdejakauma */
	      
	      for (n = 0; n < nmg*ns; n++)
		src[n] = Iin[n]/sum;
	    }
	  else if (2 == 1)
	    {
	      /* Lasketaan nettovirta */
	      
	      sum = 0;
	      for (n = 0; n < nmg*ns; n++)
		sum = sum + fabs(Iout[n] - Iin[n]);
	      
	      /* Tarkistetaan */
	      
	      if (sum < ZERO)
		Die(FUNCTION_NAME, "Jotain häikkää");
	      
	      /* Lasketaan lähdejakauma */
	      
	      for (n = 0; n < nmg*ns; n++)
		src[n] = fabs(Iout[n] - Iin[n])/sum;
	    }
	  else
	    {
	      /* Lasketaan kokonaisvuo */
	      
	      sum = 0;
	      for (n = 0; n < nmg*ns; n++)
		sum = sum + alb[n];
	      
	      /* Tarkistetaan */
	      
	      if (sum < ZERO)
		Die(FUNCTION_NAME, "Jotain häikkää");
	      
	      /* Lasketaan lähdejakauma */
	      
	      for (n = 0; n < nmg*ns; n++)
		src[n] = alb[n]/sum;
	    }

	  for (n = 0; n < nmg*ns; n++)
	    printf("%E\n", src[n]);

	  /* Populaation koko ja bätsien määrä */

	  pop = 50000;
	  nbatch = 20;

	  /* Normeerauskerroin */

	  norm = sum/(pop*nbatch*RDB[DATA_NHIST_CYCLE]);

	  /* Resetoidaan muuttujat kääntäjän hämäämiseksi */

	  x = 0.0;
	  y = 0.0;
	  z = 0.0;

	  u = 0.0;
	  v = 0.0;
	  w = 0.0;

	  ng = 0;
	  wgt = 0.0;

	  /* Silmukka bätsien yli */

	  for (m = 0; m < nbatch; m++)
	    {
	      printf("%ld/%ld\n", m + 1, nbatch);

	      /* Silmukka lähdeneutronien yli */

	      id = 0;	  
	      for (n = 0; n < pop; n++)
		{
		  /* Arvotaan parametrit lähdejakaumasta */
		  /*
		  Lahde2(&x, &y, &z, &u, &v, &w, &ng, &wgt, soos,id);
		  */


		  Lahde(&x, &y, &z, &u, &v, &w, &ng, &wgt, src, mu0, surf, ns, 
			type, params, nmg, id);


		  /* Haetaan rajapinta */
			  
		  DFPos(surf, x, y, z, &n1, &n2, &m1, &m2, &u0, &v0, &w0);

		  /* Lasketaan kosini */

		  if ((mu = fabs(u*u0 + v*v0 + w*w0)) < 0.1)
		    mu = 0.05;
		  
		  /* Skoorataan homogeeninen vuo pinnalla */

		  hom[ng + n1*nmg] += norm*wgt/mu;		      

		  /* Transport-silmukka */

		  while (TestSurface(surf, x, y, z, NO, id) == YES)
		    {
		      /* Arvotaan vapaamatka */
		      
		      l = -log(RandF(id))/tot[ng];
		      
		      /* Lasketaan etäisyys rajapintaan */
		      
		      d = SurfaceDistance(surf, params, type, ns, x, y, z, 
					  u, v, w, id);
		      
		      /* Verrataan */
		      
		      if (l < d)
			{
			  /* Siirrytään törmäyspisteeseen */
			  
			  x = x + l*u;
			  y = y + l*v;
			  z = z + l*w;
			  
			  /* Vuon törmäysestimaattori 20cm paksussa alueessa */
			  /* jaettuna 100 segmenttiin */

			  i = (long)(100.0*((x - params[0])/
					    (-20.0 - params[0])));
			  
			  if ((i > -1) && (i < 100))
			    dis[ng][i] = dis[ng][i] 
			      + norm*wgt/tot[ng]/20.0*100.0;
						  
			  /* Arvotaan vuorovaikutus */
			  
			  f = RandF(id)*tot[ng];
			  
			  /* Pelkästään sironta tarvitsee käsitellä, sillä */
			  /* kaikki muu on absorptiota */
			  
			  if (f < scatt[ng])
			    Sironta(&u, &v, &w, &ng, &wgt, nmg, scatt[ng], 
				    mubar[ng], smtx, pmtx, id);
			  else
			    break;
			}
		      else
			{
			  /* Siirrytään rajapinnalle */
			  
			  x = x + d*u;
			  y = y + d*v;
			  z = z + d*w;
			  
			  /* Haetaan rajapinta */
			  
			  DFPos(surf, x, y, z, &n1, &n2, &m1, &m2, 
				&u0, &v0, &w0);

			  /* Lasketaan kosini */

			  if ((mu = fabs(u*u0 + v*v0 + w*w0)) < 0.1)
			    mu = 0.05;

			  /* Skoorataan homogeeninen vuo pinnalla */
			  
			  hom[ng + n1*nmg] += norm*wgt/mu;		      
			  
			  /* Katkaistaan historia */
			  
			  break;
			}
		    }
		}
	    }

	  for (n = 0; n < nmg; n++)
	    fprintf(fp, "hom(%ld) = %E;\n", n + 1, hom[n]);

	  for (n = 0; n < nmg; n++)
	    fprintf(fp, "het(%ld) = %E;\n", n + 1, 
		    het[n]/RDB[DATA_NHIST_CYCLE]);

	  for (n = 0; n < nmg; n++)
	    for (i = 0; i < 100; i++)
	      fprintf(fp, "flx(%ld,%ld) = %E;\n", n + 1, i + 1, dis[n][i]);

	  /*******************************************************************/
	}

      /***********************************************************************/

      /* Vapautetaan muisti (tää myöhemmin pois) */

      Mem(MEM_FREE, scatt);
      Mem(MEM_FREE, tot);
      Mem(MEM_FREE, hom);
      Mem(MEM_FREE, src);
      Mem(MEM_FREE, alb);
      Mem(MEM_FREE, mu0);
      
      /* Seuraava universumi */

      gcu = NextItem(gcu);
      
      /************************************************************************/
    }
  
  fclose(fp);
  fprintf(out, "OK.\n\n");

#endif
}

/*****************************************************************************/

void Sironta(double *u, double *v, double *w, long *n0, double *wgt, long nmg,
	     double scatt, double mubar, const double *smtx, 
	     const double *pmtx, long id)
{
  long n1;
  double f, mu;

  /* Arvotaan osuus kokonaissirontavaikutusalasta */

  f = scatt*RandF(id);

  /* Haetaan energiaryhmä (noi ryhmät pitää kattoa) */

  for (n1 = 0; n1 < nmg; n1++)
    if ((f = f - smtx[(*n0)*nmg + n1]) < 0.0)
      break;

  /* Tarkistetaan */

  if (f > 0.0)
    Die(FUNCTION_NAME, "Sironnan arvonta epäonnistui %E %E", scatt, f);

  /* Lasketaan (n,xn)-reaktioiden painokerroin */

  f = pmtx[(*n0)*nmg + n1]/smtx[(*n0)*nmg + n1];
  CheckValue(FUNCTION_NAME, "f", "", f, 1.0, 2.0);

  /* Muokataan painoa */

  *wgt = *wgt*f;

  /* Arvotaan sirontakulma */

  mu = ArvoP1Mu(mubar, id);
  
  /* Pyöräytetään suuntakosineita */

  AziRot(mu, u, v, w, id);

  /* Asetetaan uusi energiaryhmä */
  
  *n0 = n1;
}

/*****************************************************************************/

void Lahde(double *x, double *y, double *z, double *u, double *v, double *w,
	   long *ng, double *wgt, double *src, double *mu0, long surf, long ns,
	   long type, const double *params, long nmg, long id)
{
  long m, n, i;
  double f, mu;

  /* Tarkistetaan pointteri */

  CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

  /* Arvotaan lähdevektorista pinta ja energia */

  f = RandF(id);
  for (m = 0; m < ns*nmg; m++)
    if ((f = f - src[m]) < 0.0)
      break;
  
  /* Tarkistetaan */

  if (f > 0.0)
    Die(FUNCTION_NAME, "Lähdepinnan arvonta epäonnistui");

  /* Puretaan lähdeindeksi pintaan ja energiaryhmään */
  
  i = (long)((double)m/((double)nmg));
  n = m - i*nmg;

  /* Tarkistetaan */
  
  CheckValue(FUNCTION_NAME, "i", "", i, 0, ns - 1);
  CheckValue(FUNCTION_NAME, "n", "", n, 0, nmg - 1);

  /* Arvotaan koordinaatit ja asetetaan suuntavektori pinnan normaalin */
  /* suuntaiseksi */

  switch(type)
    {
    case SURF_PX:
      {
	/* Asetetaan koordinaatit */

	*x = params[0] - EXTRAP_L;
	*y = 0.0;
	*z = 0.0;

	/* Asetetaan suuntavektori */

	*u = -1.0;
	*v = 0.0;
	*w = 0.0;

	break;
      }
    case SURF_PY:
      {
	/* Asetetaan koordinaatit */

	*x = 0.0;
	*y = params[0] - EXTRAP_L;
	*z = 0.0;

	/* Asetetaan suuntavektori */

	*u = 0.0;
	*v = -1.0;
	*w = 0.0;

	break;
      }
    case SURF_PZ:
      {
	/* Asetetaan koordinaatit */

	*x = 0.0;
	*y = 0.0;
	*z = params[0] - EXTRAP_L;

	/* Asetetaan suuntavektori */

	*u = 0.0;
	*v = 0.0;
	*w = -1.0;

	break;
      }
    default:
      {
	Die(FUNCTION_NAME, "Tätä pintaa ei oo vielä koodattu");
      }
    }
     
  /* Asetetaan energiaryhmä */

  *ng = n;
  
  /* Asetetaan paino */
  
  *wgt = 1.0;

  /* Arvotaan kulma suhteessa normaaliin (valkoinen reuna) */

  mu = sqrt(RandF(id));
  
  /* Pyöräytetäänkoordinaatteja satunnaisen kulman ympäri */
      
  AziRot(mu, u, v, w, id);
}

/*****************************************************************************/

double ArvoP1Mu(double mubar, long id)
{
  double f, mu;

  /* Arvotaan P1 sirontakulma (tää on KENO-V.a:n manuaalista sivulta 31) */
  
  f = RandF(id) - 1.0;
  
  if (fabs(mubar) < 0.0001)
    mu = f;
  else if (fabs(mubar) < 0.3333)
    mu = (sqrt(1.0 + 6.0*f*mubar + 9.0*mubar*mubar) - 1.0)/(3.0*mubar);
  else
    mu = f*(1.0 - fabs(mubar)) + mubar;

  /* Tarkistetaan */

  CheckValue(FUNCTION_NAME, "mu", "", mu, -1.0, 1.0);

  /* Palautetaan */

  return mu;
}

/*****************************************************************************/

void Lahde2(double *x, double *y, double *z, double *u, double *v, double *w,
	    long *ng, double *wgt, double **soos, long id)
{
  long n;
  double mu;

  /* Arvotaan lähdepiste */

  n = (long)(10000.0*RandF(id));

  /* Energiaryhmä ja kosini */

  mu = soos[0][n];
  *ng = (long)soos[1][n];

  /* Asetetaan koordinaatit */

  *x = -EXTRAP_L;
  *y = 0.0;
  *z = 0.0;
  
  /* Asetetaan suuntavektori */
  
  *u = -1.0;
  *v = 0.0;
  *w = 0.0;

  *wgt = 1.0;

  /* Pyöräytetäänkoordinaatteja satunnaisen kulman ympäri */
      
  AziRot(mu, u, v, w, id);
}
 
/*****************************************************************************/
 
long RegIdx(double x, double y, double z, double u, double v, double w,
	    long *cell, long id)
{
  long ptr, uni, gcu, i, ncol;
  
  /* Add to track counter */

  ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
  AddPrivateData(ptr, 1.0, id);

  /* Call geometry routine */

  *cell = WhereAmI(x, y, z, u, v, w, id);

   /* Get collision number */

  ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
  ncol = (long)GetPrivateData(ptr, id);

  /* Get universe pointer */

  if ((gcu = TestValuePair(DATA_GCU_PTR_UNI, ncol, id)) < VALID_PTR)
    return -1;

  /* Tää on tyhmä */

  uni = (long)RDB[gcu + GCU_PTR_UNIV];
  CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

  /* Get index */

  i = (long)RDB[uni + UNIVERSE_GCU_IDX];

  /* Return value */

  return i;
}

/*****************************************************************************/

void Sironta2(double *u, double *v, double *w, double *wgt, long *n0, long nmg,
	      double mubar, double *scattp, double *scattw, long id)
{
  double f, mu;
  long n;

  /* Sample new energy group */

  f = RandF(id);

  for (n = 0; n < nmg; n++)
    if ((f = f - scattp[n]) < 0.0)
      break;

  /* Check */

  if (n == nmg)
    {
      for (n = 0; n < nmg; n++)
	printf("%E\n", scattp[n]);
      Die(FUNCTION_NAME, "Unable to sample energy group");
    }

  /* Adjust weight */
  /*
  *wgt = *wgt*scattw[n];
  */
  if (scattw[n] < 1.0)
    Die(FUNCTION_NAME, "WHAT!?!?!?");

  /* Sample scattering cosine */

  mu = ArvoP1Mu(mubar, id);
  
  /* Rotate coordinates */

  AziRot(mu, u, v, w, id);

  /* Put new energy group */

  *n0 = n;
}

/*****************************************************************************/

void Fissio(double x, double y, double z, double wgt, double nubar, 
	    double *chi, long id)
{
  long nu, n, nn;
  double keff, u, v, w;

  /* Get cycle-keff */
  
  keff = RDB[DATA_CYCLE_KEFF];

  /* Adjust nubar */

  nubar = wgt*nubar/keff;

  /* Get integer part */

  nu = (long)nubar;

  /* Sample extra neutron */
  
  if (RandF(id) < nubar - (double)nu)
    nu++;

  /* Loop over source particles */

  for (n = 0; n < nu; n++)
    {
      /* Sample direction cosines */

      IsotropicDirection(&u, &v, &w, id);
      
      /* Get particle from stack */

      nn = FromStack(PARTICLE_TYPE_NEUTRON, id);

      /* Put values */
      
      WDB[nn + PARTICLE_X] = x;
      WDB[nn + PARTICLE_Y] = y;
      WDB[nn + PARTICLE_Z] = z;
      
      WDB[nn + PARTICLE_U] = u;
      WDB[nn + PARTICLE_V] = v;
      WDB[nn + PARTICLE_W] = w;
      
      WDB[nn + PARTICLE_E] = 0.0;
      WDB[nn + PARTICLE_WGT] = 1.0;
      WDB[nn + PARTICLE_FMTX_IDX] = -1.0;

      /* Put particle to bank */
      
      ToBank(nn, id);
    }
}

/*****************************************************************************/
