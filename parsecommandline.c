/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : parsecommandline.c                             */
/*                                                                           */
/* Created:       2010/11/21 (JLe)                                           */
/* Last modified: 2019/10/10 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Handles command line input                                   */
/*                                                                           */
/* Comments: - From Serpent 1.1.14                                           */
/*           - Volume checker rutiinin käsittely puuttuu                     */
/*           - Täältä kutsutaan exit(-1)                                     */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ParseCommandLine:"

/*****************************************************************************/

long ParseCommandLine(int argc, char **argv)
{
  long n, mode, np, ptr;
  char str[MAX_STR], fname[MAX_STR];
  FILE *fp;
  char tmpstr[MAX_STR];

  /* Reset file name */

  fname[0] = '\0';

  /* Avoid warning messages */

  np = 0;

  /***** Check number of arguments *******************************************/

  if ((argc < 2))
    {
      fprintf(outp, "\nUsage: %s <inputfile> [options]\n\n", argv[0]);
      fprintf(outp, "       Where <inputfile> is the file path for the main input ");
      fprintf(outp, "file and\n       the available options are:\n\n");

      fprintf(outp, "       -version           :  print version information and ");
      fprintf(outp, "exit\n");
      fprintf(outp, "       -replay            :  run simulation using random ");
      fprintf(outp, "number seed from a\n                             previous run\n");
      fprintf(outp, "       -his               :  run only burnup history in ");
      fprintf(outp, "coefficient calculation\n");
      fprintf(outp, "       -coe               :  run only restarts in ");
      fprintf(outp, "coefficient calculation\n");
      fprintf(outp, "       -plot              :  stop after geometry plot\n");
      fprintf(outp, "       -noplot            :  skip geometry plots\n");
      fprintf(outp, "       -checkvolumes <N>  :  calculate Monte Carlo estimates ");
      fprintf(outp, "for material\n                             volumes\n");
      fprintf(outp, "       -checkstl <N> <M>  :  check for holes and errors in ");
      fprintf(outp, "STL geometries\n                             by sampling ");
      fprintf(outp, "<M> directions in <N> points\n");
      fprintf(outp, "       -mpi <N>           :  run simulation in MPI mode using ");
      fprintf(outp, "<N> parallel\n                             tasks\n");
      fprintf(outp, "       -omp <M>           :  run simulation in OpenMP mode ");
      fprintf(outp, "using <M> parallel\n                             threads\n");
      fprintf(outp, "       -disperse          :  generate random particle or ");
      fprintf(outp, "pebble distribution\n                             files for ");
      fprintf(outp, "HTGR calculations\n");
      fprintf(outp, "       -rdep              :  read binary depletion file from ");
      fprintf(outp, "previous\n                             calculation and print ");
      fprintf(outp, "new output according to\n");
      fprintf(outp, "                             inventory list\n");
      fprintf(outp, "       -tracks <N>        :  draw particle tracks in ");
      fprintf(outp, "the geometry plots\n");
      fprintf(outp, "       -trackfile <N>     :  write particle tracks in ");
      fprintf(outp, "file\n");
      fprintf(outp, "       -comp <mat>        :  print pre-defined material ");
      fprintf(outp, "composition that can be\n");
      fprintf(outp, "                             copy-pasted into the input");
      fprintf(outp, "file\n");
      fprintf(outp, "       -elem <mat> <dens> :  decomposes elemental ");
      fprintf(outp, "composition into isotopic\n");
      fprintf(outp, "       -qp                :  quick plot mode (ignore overlaps)\n");
      fprintf(outp, "       -nofatal           :  ignore fatal errors\n");
      fprintf(outp, "       -mix               :  decompose mixtures in file\n");
      fprintf(outp, "       -norun             :  stop after processing\n");
      fprintf(outp, "       -matpos <coord>    :  print geometry and material data\n");
      fprintf(outp, "       -input             :  copy all inputs in a single file\n");
      fprintf(outp, "       -ip                :  launch interactive (command-line) ");
      fprintf(outp, "plotter\n");
      /*
     fprintf(outp, "       -spi <n>           :  calculate source point importances\n");
      */
      fprintf(outp, "\nSee also:\n\n");
      fprintf(outp, "http://serpent.vtt.fi/mediawiki/index.php/Installing_and_running_Serpent\n\n");
      exit(-1);
    }

  /***************************************************************************/

  /***** Check version request ***********************************************/

  if (!strcasecmp(argv[1], "-version"))
    {
      /* Print title and exit */

      PrintTitle();

      exit(-1);
    }

  /***************************************************************************/

  /***** Print standard composition ******************************************/

  if (!strcasecmp(argv[1], "-comp"))
    {
      /* Print composition and exit */

      if (argc > 3)
        StdComp(argv[2], argv[3]);
      else if (argc > 2)
        StdComp(argv[2], NULL);
      else
        fprintf(outp, "\nMissing composition name\n\n");

      exit(-1);
    }

  /***************************************************************************/

  /***** Decompose elements **************************************************/

  if (!strcasecmp(argv[1], "-elem"))
    {
      /* Print composition and exit */

      if (argc > 4)
        Element(-1, argv[2], argv[3], argv[4]);
      else if (argc > 3)
        Element(-1, argv[2], argv[3], NULL);
      else if (argc > 2)
        fprintf(outp, "\nMissing density\n\n");
      else
        fprintf(outp, "\nMissing element name\n\n");

      exit(-1);
    }

  /***************************************************************************/

  /***** Check xs test mode **************************************************/

  if (!strcasecmp(argv[1], "-testxs"))
    {
      /* Set file name */

      if (argc > 2)
        {
          /* Allocate memory for list */

          ptr = ReallocMem(DATA_ARRAY, 2);

          /* Put pointer */

          WDB[DATA_PTR_ACEDATA_FNAME_LIST] = (double)ptr;

          /* Put name */

          WDB[ptr++] = (double)PutText(argv[2]);

          /* Put null */

          WDB[ptr] = NULLPTR;

          /* Read file name */

          if (argc > 3)
            WDB[DATA_PTR_XSTEST_FNAME] = (double)PutText(argv[3]);
        }

      /* Test cross sections */

      TestXS();

      /* Exit subroutine */

      return -1;
    }

  /***************************************************************************/

  /***** Treat rest of the arguments *****************************************/

  /* Reset running mode and file pointer */

  mode = 0;
  fp = NULL;

  /* Loop over arguments */

  for (n = 1; n < argc; n++)
    {
      /* Test options */

      if (!strcmp(argv[n], "-replay"))
        mode |= MODE_REPLAY;
      else if (!strcmp(argv[n], "-his"))
        WDB[DATA_COEF_CALC_SPECIAL_MODE] = (double)SPECIAL_COEF_MODE_HIS_ONLY;
      else if (!strcmp(argv[n], "-nofatal"))
        WDB[DATA_TERMINATE_ON_DIE] = (double)NO;
      else if (!strcmp(argv[n], "-coe"))
        WDB[DATA_COEF_CALC_SPECIAL_MODE] = (double)SPECIAL_COEF_MODE_COE_ONLY;
      else if (!strcmp(argv[n], "-plot"))
        WDB[DATA_STOP_AFTER_PLOT] = (double)STOP_AFTER_PLOT_GEOM;
      else if (!strcmp(argv[n], "-noplot"))
        WDB[DATA_IGNORE_GEOM_PLOTS] = (double)YES;
      else if (!strcmp(argv[n], "-norun"))
        WDB[DATA_STOP_AFTER_PROCESSING] = (double)YES;
      else if (!strcmp(argv[n], "-mix"))
        WDB[DATA_DECOMPOSE_MIXTURES] = (double)YES;
      else if (!strcmp(argv[n], "-casematrix"))
        {
          /* Read parameters */

          if (argc - n < 4)
            {
              fprintf(errp, "\nMissing case matrix parameters.\n\n");
              exit(-1);
            }
          else
            {
              /* Read name */

              WDB[DATA_CASEMTX_RUN_PTR_NAME] = (double)PutText(argv[++n]);

              /* Read history index */

              WDB[DATA_CASEMTX_RUN_HIS_IDX] = (double)atoi(argv[++n]);

              /* Read coefficient index */

              WDB[DATA_CASEMTX_RUN_COE_IDX] = (double)atoi(argv[++n]);

              /* Set mode on */

              WDB[DATA_CASEMTX_RUN] = (double)YES;

              /* Check special options */

              if ((long)RDB[DATA_CASEMTX_RUN_COE_IDX] < 0)
                WDB[DATA_COEF_CALC_SPECIAL_MODE] =
                  (double)SPECIAL_COEF_MODE_HIS_ONLY;
              else if ((long)RDB[DATA_CASEMTX_RUN_COE_IDX] > 0)
                WDB[DATA_COEF_CALC_SPECIAL_MODE] =
                  (double)SPECIAL_COEF_MODE_COE_ONLY;
            }
        }
      else if (!strcmp(argv[n], "-spi"))
        {
          /* Source point importance calculation */

          if (argc < n + 2)
            {
              fprintf(errp, "\nNumber of source points not given.\n\n");
              exit(-1);
            }
          else if ((np = atol(argv[n + 1])) < 1)
            {
              fprintf(errp, "\nInvalid number of source points.\n\n");
              exit(-1);
            }
          else
            {
              WDB[DATA_SRC_IMP_CALC] = (double)YES;
              WDB[DATA_SRC_IMP_NPS] = (double)np;
            }

          n++;
        }
      else if (!strcmp(argv[n], "-tracks") || !strcmp(argv[n], "-trackfile"))
        {
          /* Track plotter */

          if (argc < n + 2)
            {
              fprintf(errp, "\nNumber of histories not given.\n\n");
              exit(-1);
            }
          else if ((np = atol(argv[n + 1])) < 1)
            {
              fprintf(errp, "\nInvalid number of histories.\n\n");
              exit(-1);
            }
          else
            {
              WDB[DATA_STOP_AT_BOUNDARY] = (double)YES;
              WDB[DATA_STOP_AFTER_PLOT] = (double)STOP_AFTER_PLOT_TRACKS;
              WDB[DATA_TRACK_PLOTTER_HIS] = (double)np;
              WDB[DATA_QUICK_PLOT_MODE] = (double)YES;
            }

          if (!strcmp(argv[n], "-trackfile"))
            WDB[DATA_TRACK_PLOTTER_FILE] = (double)YES;
          n++;
        }
      else if (!strcmp(argv[n], "-testgeom"))
        {
          fprintf(errp, "Geometry tester is not available in Serpent 2.\n");
          exit(-1);
        }
      else if (!strcasecmp(argv[n], "-rdep"))
        WDB[DATA_PARTICLE_REDEPLETE_MODE] = YES;
      else if (!strcmp(argv[n], "-disperse"))
        WDB[DATA_PARTICLE_DISPERSER_MODE] = YES;
      else if (!strcmp(argv[n], "-qp"))
        WDB[DATA_QUICK_PLOT_MODE] = YES;
      else if (!strcmp(argv[n], "-ip"))
        WDB[DATA_INTERACTIVE_PLOT_MODE] = YES;
      else if (!strcmp(argv[n], "-input"))
        WDB[DATA_COPY_INPUTS] = YES;
      else if (!strcmp(argv[n], "-matpos"))
        {
          /* Check if input file name is given */

          if (fname[0] == '\0')
            {
              /* Option before input file, read file name */

              WDB[DATA_MATPOS_PTR_FILE] = (double)PutText(argv[++n]);

              /* Set number of points */

              WDB[DATA_MATPOS_N_PTS] = -1.0;
            }
          else
            {
              /* Assume option is given last */

              if ((np  = argc - n - 1) % 3)
                Error(0, "Invalid number of coordinates given");

              /* Update index */

              n++;

              /* Put number of points */

              WDB[DATA_MATPOS_N_PTS] = (double)(np/3);

              /* Allocate memory for points */

              ptr = ReallocMem(DATA_ARRAY, np);
              WDB[DATA_MATPOS_PTR_COORD] = (double)ptr;

              /* Loop over data */

              while (n < argc)
                WDB[ptr++] = atof(argv[n++]);
            }
        }
      else if ((!strcmp(argv[n], "-checkvolume")) ||
               (!strcmp(argv[n], "-checkvolumes")))
        {
          /* Material volume test. */

          if (argc < n + 2)
            {
              fprintf(errp, "Number of test points not given.\n");
              exit(-1);
            }
          else
            {
              /* Get number of random points */

              WDB[DATA_VOLUME_MC_NMAX] = atof(argv[++n]);

              /* Set running mode */

              WDB[DATA_VOLUME_CALCULATION_MODE] = YES;
              WDB[DATA_VOLUME_MC] = YES;
            }
        }
      else if (!strcmp(argv[n], "-checkstl"))
        {
          /* STL geometry test */

          if (argc < n + 3)
            {
              fprintf(errp, "Number of test points not given.\n");
              exit(-1);
            }
          else
            {
              /* Get number of random points and directions */

              WDB[DATA_STL_TEST_N_PTS] = atof(argv[++n]);
              WDB[DATA_STL_TEST_N_DIR] = atof(argv[++n]);
            }
        }

#ifdef MPI

      else if (!strcmp(argv[n], "-mpi"))
        {
          mode |= MODE_MPI;

          /* Get number of tasks */

          if (argc < n + 2)
            {
              fprintf(errp, "Number of MPI tasks not given.\n");
              exit(-1);
            }
          else if ((np = atoi(argv[n + 1])) < 2)
            {
              fprintf(errp, "Invalid MPI task number \"%s\".\n",argv[n + 1]);
              exit(-1);
            }
          n++;
        }

#else

      else if (!strcmp(argv[n], "-mpi"))
        {
          fprintf(errp,
                  "\nMPI parallel calculation not available in compilation.\n\n");
          exit(-1);
        }

#endif

#ifdef OPEN_MP

      else if (!strcmp(argv[n], "-omp"))
        {
          /* Get number of threads */

          if (argc < n + 2)
            {
              fprintf(errp, "Number of OpenMP threads not given.\n");
              exit(-1);
            }
          else if (!strcmp(argv[n + 1], "max"))
            {
              /* Use maximum available number of threads */

              WDB[DATA_OMP_MAX_THREADS] = -1.0;
            }
          else if (atoi(argv[n + 1]) < 1)
            {
              fprintf(errp, "Invalid OpenMP thread number \"%s\".\n",
                      argv[n + 1]);
              exit(-1);
            }
          else if (atoi(argv[n + 1]) > MAX_OMP_THREADS)
            {
              fprintf(errp, "Maximum number of OpenMP threads is %d",
                      MAX_OMP_THREADS);
              exit(-1);
            }
          else
            {
              /* Set number of threads */

              WDB[DATA_OMP_MAX_THREADS] = (double)atoi(argv[n + 1]);
            }

          n++;
        }

#else

      else if (!strcmp(argv[n], "-omp"))
        {
          fprintf(errp, "OpenMP mode not available in compilation.\n");
          exit(-1);
        }

#endif

      else if (fp == NULL)
        {

          /* Argument is a file name, exit if disperse mode */

          if ((long)RDB[DATA_PARTICLE_DISPERSER_MODE] == YES)
            {
              WDB[DATA_PTR_INPUT_FNAME] = (double)PutText(argv[n]);
              return 0;
            }

          /* Check that file exists */

          else if ((fp = fopen(argv[n], "r")) == NULL)
            {
              fprintf(errp, "Input file \"%s\" does not exist.\n", argv[n]);
              exit(-1);
            }

          /* Remember file name */

          strcpy(fname, argv[n]);
        }
      else
        {
          /* One file name already read */

          fprintf(errp, "Invalid parameter \"%s\".\n", argv[n]);
          exit(-1);
        }
    }

  /***************************************************************************/

  /***** Check that input file is given **************************************/

  /* Exit if disperse mode */

  if ((long)RDB[DATA_PARTICLE_DISPERSER_MODE] == YES)
    return 0;

  if (fp == NULL)
    {
      fprintf(errp, "No input file given.\n");
      exit(-1);
    }
  else
    fclose(fp);

  /***************************************************************************/

  /***** Check for MPI-mode **************************************************/

#ifdef MPI

  if (mode & MODE_MPI)
    {
      /* Reset string */

      tmpstr[0] = '\0';

      /* Parse command string */

      sprintf(str, "%s -np %ld %s %s ", MPIRUN_PATH, np, argv[0], fname);

      /* Add some special options */

      if (mode & MODE_REPLAY)
        strcat(str, "-replay ");

      if ((long)RDB[DATA_OMP_MAX_THREADS] == -1)
        strcat(str, "-omp max ");
      else if ((long)RDB[DATA_OMP_MAX_THREADS] > 1)
        {
          sprintf(tmpstr, "-omp %ld ", (long)RDB[DATA_OMP_MAX_THREADS]);
          strcat(str, tmpstr);
        }

      /* Call recursively */

      if (system(str) == -1)
        Die(FUNCTION_NAME, "Recursive call failed");

      /* Free memory */

      FreeMem();

      /* Finalize MPI */

      FinalizeMPI();

      /* Exit */

      return -1;
    }

#endif

  /***************************************************************************/

  /***** Check that seed file exists if in replay mode ***********************/

  if (mode & MODE_REPLAY)
    {
      sprintf(str, "%s.seed", fname);

      if ((fp = fopen(str, "r")) != NULL)
        {
          if (fscanf(fp, "%lu", &parent_seed) == EOF)
            Die(FUNCTION_NAME, "fscanf error");

          fclose(fp);
        }
      else
        {
          fprintf(errp, "Seed file \"%s\" does not exist.\n", str);
          exit(-1);
        }

      /* Set replay option */

      WDB[DATA_OPTI_REPLAY] = (double)YES;
    }
  else
    {
      /* Write seed file */

      sprintf(str, "%s.seed", fname);

      if ((fp = fopen(str, "w")) == NULL)
        {
          fprintf(errp, "%s Unable to open seed file for writing.\n\n",
                  FUNCTION_NAME);
          exit(-1);
        }

      fprintf(fp, "%lu\n", parent_seed);

      fclose(fp);
    }

  /***************************************************************************/

  /* Set file name */

  WDB[DATA_PTR_INPUT_FNAME] = (double)PutText(fname);

  /* Remove open output file */

  if ((long)RDB[DATA_COPY_INPUTS] == YES)
    {
      sprintf(tmpstr, "%s.input", fname);
      remove(tmpstr);
    }

  /* Init random number generator */

  srand48(parent_seed);

  /* Exit subroutine */

  return 0;
}

/*****************************************************************************/
