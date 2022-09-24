/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : testxs.c                                       */
/*                                                                           */
/* Created:       2010/12/26 (JLe)                                           */
/* Last modified: 2018/02/02 (JLe)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Test routine for reading and processing cross sections etc.  */
/*                                                                           */
/* Comments: - Toimii vain yhdellä xsdata directory filellä                  */
/*           - En testannut, toimiiko (TVi 2015-05-19)                       */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "TestXS:"

#define TOFILE

void TestDistribution(long, double, long);
void TestDistributions(long, double, long);

/*****************************************************************************/

void TestXS()
{
  long nmat, nnuc, nace, n, m, i, ace, ptr, mat, iso, count, N, poi, ii;
  long loc, mix, nth, tl, th, tg, id, nucmin, nucmax, matmin, matmax, fiss;
  char tmpstr[MAX_STR], xsdata[MAX_STR], **nuc, fname[MAX_STR];
  double E;
  FILE *inpfile;
  unsigned long seed;

  fprintf(outp, "\nRunning cross section tester routine...\n\n");

  /* Get file name */

  if ((long)RDB[DATA_PTR_XSTEST_FNAME] > VALID_PTR)
    sprintf(fname, "%s", GetText(DATA_PTR_XSTEST_FNAME));
  else
    sprintf(fname, "tester");

  /* Number of test points */

  N = 10000;

  /* Set maximum materials and nuclides per material */

  nucmin = 2;
  nucmax = 50;
  matmin = 2;
  matmax = 5;

  /***************************************************************************/

  /***** Read available nuclide names ****************************************/

  /* Read ace directory file */

  ReadDirectoryFile();

  /* Loop over ace files and count */

  nace = 0;
  ace = (long)RDB[DATA_PTR_ACE0];

  while (ace > VALID_PTR)
    {
      /* Add counter */

      nace++;

      /* Next */

      ace = (long)ACE[ace + ACE_PTR_NEXT];
    }

  /* Allocate memory */

  nuc = (char **)Mem(MEM_ALLOC, nace, sizeof(char *));

  for(n = 0; n < nace; n++)
    nuc[n] = (char *)Mem(MEM_ALLOC, 200, sizeof(char));

  /* Read data */

  n = 0;
  ace = (long)RDB[DATA_PTR_ACE0];

  while (ace > VALID_PTR)
    {
      /* Copy name */

      if((ACE[ace + ACE_TYPE] == NUCLIDE_TYPE_SAB) ||
         (ACE[ace + ACE_TEMP] == 300.0))
        {
          WDB[DATA_DUMMY] = ACE[ace + ACE_PTR_NAME];
          sprintf(nuc[n++], "%s", GetText(DATA_DUMMY));
        }

      /* Next */

      ace = (long)ACE[ace + ACE_PTR_NEXT];
    }

  nace = n;

  /* Remember directory file name */

  ptr = (long)RDB[DATA_PTR_ACEDATA_FNAME_LIST];
  sprintf(xsdata, "%s", GetText(ptr));

  /***************************************************************************/

  /***** Loop ****************************************************************/

  /* Loop over repeats and counts */

  for (count = 0; count < 100000; count++)
    {
      /***********************************************************************/

      /***** Generate input file for calculations ****************************/

      /* Remove old output file */

      sprintf(tmpstr, "%s.out", fname);
      remove(tmpstr);

      /* Free memory */

      FreeMem();

      /* Init data array */

      InitData();

      /* Reset thermal library counter */

      nth = 0;

      /* Reset fissile flag */

      fiss = NO;

      /* Allocate memory for directory file list */

      ptr = ReallocMem(DATA_ARRAY, 2);

      /* Put pointer */

      WDB[DATA_PTR_ACEDATA_FNAME_LIST] = (double)ptr;

      /* Put name */

      WDB[ptr++] = (double)PutText(xsdata);

      /* Put null */

      WDB[ptr] = NULLPTR;

      /* Init random number sequence */

      seed = ReInitRNG(id);
      SEED[id*RNG_SZ] = seed;

      id = 0;

      /* Sample number of materials */

      nmat = (long)(RandF(id)*(matmax - matmin) + matmin);

      /* Loop over materials */

      for (n = 0; n < nmat; n++)
        {
          /* Create material */

          mat = NewItem(DATA_PTR_M0, MATERIAL_BLOCK_SIZE);

          /* Set used-flag */

          SetOption(mat + MATERIAL_OPTIONS, OPT_USED);

          /* Put name */

          sprintf(tmpstr, "mat%ld", n + 1);
          WDB[mat + MATERIAL_PTR_NAME] = (double)PutText(tmpstr);

          /* Density */

          WDB[mat + MATERIAL_ADENS] = RandF(id);

          /* Reset thermal flags */

          tl = 0;
          th = 0;
          tg = 0;

          /* Sample number of nuclides per material */

          nnuc = (long)(RandF(id)*(nucmax - nucmin) + nucmin);

          /* Loop over composition */

          for (m = 0; m < nnuc; m++)
            {
              /* Create new nuclide*/

              iso = NewItem(mat + MATERIAL_PTR_COMP, COMPOSITION_BLOCK_SIZE);

              /* Sample nuclide */

              do
                {
                  i = (long)(RandF(id)*nace);

                  /* Test thermal scattering */

                  if ((nuc[i][0] == 'l') && (tl == 0))
                    {
                      /* Put flag */

                      tl++;

                      /* Put nuclide name */

                      WDB[iso + COMPOSITION_PTR_NUCLIDE] =
                        (double)PutText("1001.03c");

                      /* Sample fraction */

                      WDB[iso + COMPOSITION_ADENS] = RandF(id);

                      /* S(a,b) name */

                      sprintf(tmpstr, "lwtr%ld", nth++);

                      /* Add S(a,b) data */

                      ptr = NewItem(DATA_PTR_T0, THERM_BLOCK_SIZE);
                      WDB[ptr + THERM_PTR_ALIAS] =
                        (double)PutText(tmpstr);

                      loc = NewItem(ptr + THERM_PTR_SAB, SAB_BLOCK_SIZE);
                      WDB[loc + SAB_PTR_NAME] = (double)PutText(nuc[i]);

                      WDB[ptr + THERM_T] = -1.0;

                      /* Add material entry */

                      ptr = NewItem(mat + MATERIAL_PTR_SAB,
                                    THERM_BLOCK_SIZE);
                      WDB[ptr + THERM_PTR_ALIAS] =
                        (double)PutText(tmpstr);
                      WDB[ptr + THERM_ZA] = 1001.0;
                    }
                  else if ((nuc[i][0] == 'h') && (th == 0))
                    {
                      /* Put flag */

                      th++;

                      /* Put nuclide name */

                      WDB[iso + COMPOSITION_PTR_NUCLIDE] =
                        (double)PutText("1002.03c");

                      /* Sample fraction */

                      WDB[iso + COMPOSITION_ADENS] = RandF(id);

                      /* S(a,b) name */

                      sprintf(tmpstr, "hwtr%ld", nth++);

                      /* Add S(a,b) data */

                      ptr = NewItem(DATA_PTR_T0, THERM_BLOCK_SIZE);
                      WDB[ptr + THERM_PTR_ALIAS] =
                        (double)PutText(tmpstr);

                      loc = NewItem(ptr + THERM_PTR_SAB, SAB_BLOCK_SIZE);
                      WDB[loc + SAB_PTR_NAME] = (double)PutText(nuc[i]);

                      WDB[ptr + THERM_T] = -1.0;

                      /* Add material entry */

                      ptr = NewItem(mat + MATERIAL_PTR_SAB,
                                    THERM_BLOCK_SIZE);
                      WDB[ptr + THERM_PTR_ALIAS] =
                        (double)PutText(tmpstr);
                      WDB[ptr + THERM_ZA] = 1002.0;
                    }
                  else if ((nuc[i][0] == 'g') && (tg == 0))
                    {
                      /* Put flag */

                      tg++;

                      /* Put nuclide name */

                      WDB[iso + COMPOSITION_PTR_NUCLIDE] =
                        (double)PutText("6000.03c");

                      /* Sample fraction */

                      WDB[iso + COMPOSITION_ADENS] = RandF(id);

                      /* S(a,b) name */

                      sprintf(tmpstr, "grph%ld", nth++);

                      /* Add S(a,b) data */

                      ptr = NewItem(DATA_PTR_T0, THERM_BLOCK_SIZE);
                      WDB[ptr + THERM_PTR_ALIAS] =
                        (double)PutText(tmpstr);

                      loc = NewItem(ptr + THERM_PTR_SAB, SAB_BLOCK_SIZE);
                      WDB[loc + SAB_PTR_NAME] = (double)PutText(nuc[i]);

                      WDB[ptr + THERM_T] = -1.0;

                      /* Add material entry */

                      ptr = NewItem(mat + MATERIAL_PTR_SAB,
                                    THERM_BLOCK_SIZE);
                      WDB[ptr + THERM_PTR_ALIAS] =
                        (double)PutText(tmpstr);
                      WDB[ptr + THERM_ZA] = 6000.0;
                    }
                  else if (nuc[i][strlen(nuc[i]) - 1] == 'c')
                    {
                      /* Re-sample id */

                      if ((ii = RandF(0)*7) == 0)
                        sprintf(&nuc[i][strlen(nuc[i]) - 3], "00c");
                      else if (ii == 1)
                        sprintf(&nuc[i][strlen(nuc[i]) - 3], "03c");
                      else if (ii == 2)
                        sprintf(&nuc[i][strlen(nuc[i]) - 3], "06c");
                      else if (ii == 3)
                        sprintf(&nuc[i][strlen(nuc[i]) - 3], "09c");
                      else if (ii == 4)
                        sprintf(&nuc[i][strlen(nuc[i]) - 3], "12c");
                      else if (ii == 5)
                        sprintf(&nuc[i][strlen(nuc[i]) - 3], "15c");
                      else if (ii == 6)
                        sprintf(&nuc[i][strlen(nuc[i]) - 3], "18c");
                      else
                        Die(FUNCTION_NAME, "Error");

                      /* Put nuclide name */

                      WDB[iso + COMPOSITION_PTR_NUCLIDE] =
                        (double)PutText(nuc[i]);

                      /* Sample fraction */

                      WDB[iso + COMPOSITION_ADENS] = RandF(id);
                    }
                  else
                    i = -1;
                }
              while (i == -1);

              /* Check fissile (check U and Pu only) */

              if (nuc[i][0] == '9')
                if ((nuc[i][1] == '2') || (nuc[i][1] == '4'))
                  fiss = YES;
            }
        }

      /* Create mixture */

      mat = NewItem(DATA_PTR_M0, MATERIAL_BLOCK_SIZE);

      /* Set used-flag */

      SetOption(mat + MATERIAL_OPTIONS, OPT_USED);

      /* Put name */

      sprintf(tmpstr, "mixer");
      WDB[mat + MATERIAL_PTR_NAME] = (double)PutText(tmpstr);

      /* Put first material */

      ptr = (long)RDB[DATA_PTR_M0];

      mix = NewItem(mat + MATERIAL_PTR_MIX, MIXTURE_BLOCK_SIZE);
      WDB[mix + MIXTURE_PTR_MAT] = RDB[ptr + MATERIAL_PTR_NAME];
      WDB[mix + MIXTURE_VFRAC] = RandF(id);

      /* Loop over others */

      n = 0;

      do
        {
          ptr = (long)RDB[DATA_PTR_M0];
          while (ptr < mat)
            {
              if (RandF(id) < 0.5)
                {
                  mix = NewItem(mat + MATERIAL_PTR_MIX,
                                MIXTURE_BLOCK_SIZE);
                  WDB[mix + MIXTURE_PTR_MAT] =
                    RDB[ptr + MATERIAL_PTR_NAME];
                  WDB[mix + MIXTURE_VFRAC] = RandF(id);

                  n++;
                }

              ptr = NextItem(ptr);
            }
        }
      while (n == 0);

      /* Open file for writing */

      inpfile = fopen(fname, "w");

      /* Check fissile and add equilibrium xenon / samarium calculation */
      /* 2.11.2015: Toi ei toimi jostain syystä */

      poi = NO;

      if (fiss == YES)
        {
          if (RandF(0) < -0.1)
            {
              fprintf(inpfile, "set xenon 1\n");
              poi = YES;
            }

          if (RandF(0) < -0.1)
            {
              fprintf(inpfile, "set samarium 1\n");
              poi = YES;
            }
        }

      /* Create geometry */

      fprintf(inpfile, "pin 0\n");

      for (n = 0; n < nmat; n++)
        fprintf(inpfile, "mat%ld %ld\n", n + 1, n + 1);

      if (poi == NO)
        fprintf(inpfile, "mixer %ld\n", n + 1);

      fprintf(inpfile, "outside\n");

      /* Create source */

      fprintf(inpfile, "set nps 10000\n");
      fprintf(inpfile, "src 1 sp 0 0 0\n");
      fprintf(inpfile, "set gcut 5\n");

      /* Loop over materials */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Check if material is material or mixture */

          if ((iso = (long)RDB[mat + MATERIAL_PTR_COMP]) > VALID_PTR)
            {
              fprintf(inpfile, "\nmat %s %E ",
                      GetText(mat + MATERIAL_PTR_NAME),
                      RDB[mat + MATERIAL_ADENS]);

              ptr = (long)RDB[mat + MATERIAL_PTR_SAB];
              while (ptr > VALID_PTR)
                {
                  fprintf(inpfile, "moder %s %ld ",
                          GetText(ptr + THERM_PTR_ALIAS),
                          (long)RDB[ptr + THERM_ZA]);

                  ptr = NextItem(ptr);
                }

              if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
                fprintf(inpfile, "burn 1\n");
              else
                fprintf(inpfile, "\n");

              /* Loop over composition */

              while(iso > VALID_PTR)
                {
                  if (RDB[iso + COMPOSITION_ADENS] > 100.0)

                    fprintf(inpfile, "%s %ld\n",
                            GetText(iso + COMPOSITION_PTR_NUCLIDE),
                            (long)RDB[iso + COMPOSITION_ADENS]);
                  else
                    fprintf(inpfile, "%s %f\n",
                            GetText(iso + COMPOSITION_PTR_NUCLIDE),
                            RDB[iso + COMPOSITION_ADENS]);

                  /* Next isotope */

                  iso = NextItem(iso);
                }
            }
          else if (((mix = (long)RDB[mat + MATERIAL_PTR_MIX]) > VALID_PTR) &&
                   (poi == NO))
            {
              fprintf(inpfile, "\nmix %s\n",
                      GetText(mat + MATERIAL_PTR_NAME));

              /* Loop over composition */

              while(mix > VALID_PTR)
                {
                  fprintf(inpfile, "%s %E\n",
                          GetText(mix + MIXTURE_PTR_MAT),
                          RDB[mix + MIXTURE_VFRAC]);

                  /* Next isotope */

                  mix = NextItem(mix);
                }
            }

          /* Next material */

          mat = NextItem(mat);
        }

      fprintf(inpfile, "\n");

      ptr = (long)RDB[DATA_PTR_T0];
      while (ptr > VALID_PTR)
        {
          loc = (long)RDB[ptr + THERM_PTR_SAB];

          fprintf(inpfile, "therm %s %s\n",
                  GetText(ptr + THERM_PTR_ALIAS),
                  GetText(loc + SAB_PTR_NAME));

          ptr = NextItem(ptr);
        }

      /* Libraries */

      fprintf(inpfile, "\nset acelib \"%s\"\n", xsdata);
      fprintf(inpfile, "\nset declib \"sss_jeff31.dec\"\n");
      fprintf(inpfile, "\nset nfylib \"sss_jeff31.nfy\"\n");

      if (RandF(0) < 0.2)
        fprintf(inpfile, "\nset egrid %E\n", pow(10.0, -(RandF(0)*3.0 + 2.0)));

      /*
      if (RDB[DATA_ERG_TOL] > 0.0)
        fprintf(inpfile, "\nset egrid %E %E %E\n", RDB[DATA_ERG_TOL],
                -RDB[DATA_NEUTRON_EMIN], -RDB[DATA_NEUTRON_EMAX]);
      */

      fprintf(inpfile, "\nset opti %ld\n", (long)(4.0*RandF(id)) + 1);
      fprintf(inpfile, "set ures %ld\n", (long)(2.0*RandF(id)));

      if (RandF(0) < 0.5)
        fprintf(inpfile, "set dix 1\n");

      fprintf(inpfile, "set seed %lu\n", parent_seed);

      /*
      fprintf(inpfile, "\nset xsplot 1000 %E %E\n", -RDB[DATA_NEUTRON_EMIN],
              -RDB[DATA_NEUTRON_EMAX]);
      */
      fclose(inpfile);

      /***********************************************************************/

      /***** Read data, process and test *************************************/

      /* Free memory */

      FreeMem();

      /* Init data array */

      InitData();

#ifdef TOFILE

      /* Open output file */

      sprintf(tmpstr, "%s.log", fname);
      outp = fopen(tmpstr, "w");

      /* Set line-buffering for stdout */

      setlinebuf(outp);

#endif

      errp = outp;

      /* Init OpenMP related stuff */

      InitOMP();

      /* Put input file name */

      WDB[DATA_PTR_INPUT_FNAME] = (double)PutText(fname);

      /* Read input */

      ReadInput(GetText(DATA_PTR_INPUT_FNAME));
      fprintf(outp, "\n");

      /* Set optimization */

      SetOptimization();

      /* Set some pointers and modes */

      WDB[DATA_MICRO_PTR_EGRID] = -1.0;
      WDB[DATA_OPTI_GC_CALC] = (double)NO;

      /* Set material-wise flags */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Set flag */

          SetOption(mat + MATERIAL_OPTIONS, OPT_USED);
          SetOption(mat + MATERIAL_OPTIONS, OPT_INCLUDE_MAJORANT);

          /* Next */

          mat = NextItem(mat);
        }

      /* Set material pointers */

      FindMaterialPointers();

      /* Process nuclides */

      ProcessNuclides();

      /* Process energy grids */

      UnionizeGrid();

      /* Process XS data */

      ProcessXSData();

      /* Generate cache-optimized block */

      CacheXS();

      /* Process materials */

      ProcessMaterials();

      /* Expand PRIVA, BUF and RES2 arrays for OpenMP parallel calculation */

      ExpandPrivateArrays();

      /* Disallow memory allocation from here on*/

      Mem(MEM_DENY);

      /* Prepare transport cycle */

      PrepareTransportCycle();

      fprintf(outp, "\nLoop %ld, optimization mode %ld, mem %1.2f Mb...\n\n",
              count + 1, (long)RDB[DATA_OPTI_MODE],
              RDB[DATA_TOTAL_BYTES]/MEGA);

      fprintf(outp, "Testing with %ld samples...\n", N);

      /* Reset OpenMP thread number */

      id = 0;

      /* Loop over samples */

      for (n = 0; n < N; n++)
        {
          /* Init random number sequence */

          seed = ReInitRNG(n);
          SEED[id*RNG_SZ] = seed;

          /* Sample target material */

          m = (long)((double)nmat*RandF(id));

          mat = (long)RDB[DATA_PTR_M0];
          for (i = 0; i < m; i++)
            mat = NextItem(mat);

          /* Check pointer */

          CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

          /* Sample energy */

          E = exp(RandF(id)*(log(RDB[DATA_NEUTRON_EMAX]) -
                             log(RDB[DATA_NEUTRON_EMIN]))
                  + log(RDB[DATA_NEUTRON_EMIN]));

          /* Reaction list sums */

          CheckReaListSum(mat, PARTICLE_TYPE_NEUTRON, E, NO, id);

          /* Sample reaction */

          ptr = SampleReaction(mat, PARTICLE_TYPE_NEUTRON, E, 1.0, id);

          /* Test distribution (ensimmäinen tapa sämplää reaktiot  */
          /* niiden todennäköisyyden mukaan, jälkimmäinen käy läpi */
          /* kaikki) */

          if (1 == 2)
            TestDistribution(ptr, E, id);
          else
            {
              /* Loop over composition */
              
              iso = (long)RDB[mat + MATERIAL_PTR_COMP];
              while (iso > VALID_PTR)
                {
                  /* Pointer to nuclide */
                  
                  ptr = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                  
                  /* Test distributions */
                  
                  TestDistributions(ptr, E, id);
                  
                  /* Next nuclide */
                  
                  iso = NextItem(iso);
                }
            }
        }

      fprintf(outp, "OK.\n\n");

#ifdef TOFILE

      /* Close file */

      fclose(outp);

#endif

      /*******************************************************************/
    }

   /*************************************************************************/

  /* Free memory */

  FreeMem();

  /* Free nuclide array */

  for(n = 0; n < nace; n++)
    Mem(MEM_FREE, nuc[n]);

  Mem(MEM_FREE, nuc);
}

/*****************************************************************************/

void TestDistribution(long rea, double E, long id)
{
  long mt;
  double Eout, mu, u, v, w;

  /* Check energy boundaries */

  if ((E < RDB[rea + REACTION_EMIN]) || (E > RDB[rea + REACTION_EMAX]))
    return;
  
  /* Sample random value for mu to prevent warnings */
  
  mu = RandF(id);
  IsotropicDirection(&u, &v, &w, id);
  
  /* Get mt */
  
  mt = (long)RDB[rea + REACTION_MT];
  
  /* Check */
  
  if ((mt == 1002) || (mt == 1004))
    {
      /* Copy energy */
      
      Eout = E;
      
      /* Sample S(a,b) scattering */
      
      SabScattering(rea, &Eout, &u, &v, &w, id);
    }
  else if ((long)RDB[rea + REACTION_PTR_ERG] > VALID_PTR)
    {
      /* Sample energy distribution */
      
      SampleENDFLaw(rea, -1, E, &Eout, &mu, id);
    }
  
  /* Sample angular distribution */
  
  mu = SampleMu(rea, -1, E, NULL, NULL, id);
  CheckValue(FUNCTION_NAME, "mu", "", mu, -1.0, 1.0);
}

/*****************************************************************************/

void TestDistributions(long nuc, double E, long id)
{
  long rls, rea, mt;
  double Eout, mu, u, v, w;

  /* Loop over reactions */

  rls = (long)RDB[nuc + NUCLIDE_PTR_SAMPLE_REA_LIST];
  while (rls > VALID_PTR)
    {
      /* Pointer to reaction data */

      rea = (long)RDB[rls + RLS_DATA_PTR_REA];
      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

      /* Sample energy between boundaries */

      E = exp(RandF(id)*(log(RDB[rea + REACTION_EMAX]) -
                         log(RDB[rea + REACTION_EMIN]))
              + log(RDB[rea + REACTION_EMIN]));

      /* Sample random value for mu to prevent warnings */

      mu = RandF(id);
      IsotropicDirection(&u, &v, &w, id);

      /* Get mt */

      mt = (long)RDB[rea + REACTION_MT];

      /* Check */

      if ((mt == 1002) || (mt == 1004))
        {
          /* Copy energy */

          Eout = E;

          /* Sample S(a,b) scattering */

          SabScattering(rea, &Eout, &u, &v, &w, id);
        }
      else if ((long)RDB[rea + REACTION_PTR_ERG] > VALID_PTR)
        {
          /* Sample energy distribution */

          SampleENDFLaw(rea, -1, E, &Eout, &mu, id);
        }

      /* Sample angular distribution */

      mu = SampleMu(rea, -1, E, NULL, NULL, id);
      CheckValue(FUNCTION_NAME, "mu", "", mu, -1.0, 1.0);

      /* Next reaction */

      rls = NextItem(rls);
    }
}

/*****************************************************************************/
