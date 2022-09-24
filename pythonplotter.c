/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : pythonplotter.c                                */
/*                                                                           */
/* Created:       2019/10/09 (JLe)                                           */
/* Last modified: 2020/06/15 (SÄk)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Interface routines for interactive python plotter.           */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#define PYTHONPLOTTER_VERSION 11


#include "header.h"
#include "locations.h"
#include <stdio.h>

#define FUNCTION_NAME "PythonPlotter:"

/* This struct conveys the geometry limits */

struct geom_stats
{
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double zmin;
  double zmax;
};

/* This struct conveys a material */

struct material
{
  double density; /* g/cm3 */
  long color[3]; /* R,G,B */
  long numberOfNuclides; /* For a future implementation of a routine that reads the nuclide list.*/
  char name[MAX_STR]; /* Null terminated */
};

/* This struct is used convey information about a x,y,z point from serpent MatPos to python */

struct pinfo {
  double x;
  double y;
  double z;
  double time;
  long   universe;
  char   universe_name[MAX_STR];
  long   cell;
  char   cell_name[MAX_STR];
  long   material;
  char   material_name[MAX_STR];
  long   material_color_RGB[3];
};
typedef struct pinfo pinfo;


long getVersion(char *serpentVersion, long *pythonplotterversion){

       char  *serpver =CODE_VERSION;

       *pythonplotterversion = PYTHONPLOTTER_VERSION;


       strncpy(serpentVersion,serpver,MAX_STR);

       return 0;

}



#define RGB(R,G,B) R*1000000 + G* 1000 + B

/*****************************************************************************/

long MatPosUCM(double x, double y, double z, double t,
               long *universe, long *cell, long *material,
               char **universe_name, char **cell_name, char **material_name,
               long *mat_color_R, long *mat_color_G, long *mat_color_B )
{

  double u,v,w;
  double f,T;
  long rgb,r,g,b;

  long rgb_void = 0;
  long rgb_nocell = RGB(0,255,0);

  /*
   *
   long rgb_mucell = RGB(255,0,0);
   long rgb_perror = RGB(255,255,0);
   long rgb_denser = RGB(255,0,255);
  */

  /* Sample direction (this is necessary for STL geometries) */

  IsotropicDirection(&u, &v, &w, 0);

  /* Find position */

  if ((*cell = WhereAmI(x, y, z, u, v, w, 0)) < 0)
    {
      *universe = -1;
      *universe_name = "N/A";
      *cell_name = "Geometry error";
      *material = -1;
      *material_name = "N/A";
      rgb = rgb_nocell;
      r = (long)(rgb/1000000.0);
      g = (long)((rgb - 1000000.0*r)/1000.0);
      b = (long)(rgb - 1000000.0*r - 1000.0*g);

      *mat_color_R = r;
      *mat_color_G = g;
      *mat_color_B = b;

      return 0;
    }

  /* Get universe */

  *universe = (long)RDB[*cell + CELL_PTR_UNI];
  CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, *universe);

  /* Reset density and temperature */

  f = 1.0;
  T = 0.0;

  /* Check if cell has material */

  if ((*material = (long)RDB[*cell + CELL_PTR_MAT]) > VALID_PTR)
    {
      /* Get material pointer */

      *material = MatPtr(*material, 0);

      /* Get point from interface */

      IFCPoint(*material, &f, &T, -1.0, 0);
    }

  /* Print or save to the struct */

  if ((long)RDB[*universe + UNIVERSE_PTR_NAME] < VALID_PTR)
    *universe_name = "N/A";
  else
    *universe_name = GetText(*universe + UNIVERSE_PTR_NAME);

  *cell_name = GetText(*cell + CELL_PTR_NAME);

  if ((long)RDB[*cell + CELL_TYPE] == CELL_TYPE_OUTSIDE)
    {
      *material_name = "outside";
      rgb = rgb_void;
    }
  else if (*material < VALID_PTR)
    {
      *material_name = "void";
      rgb = rgb_void;
    }
  else
    {
      *material_name =  GetText(*material + MATERIAL_PTR_NAME);
      rgb = (long)RDB[*material + MATERIAL_RGB];
    }

  r = (long)(rgb/1000000.0);
  g = (long)((rgb - 1000000.0*r)/1000.0);
  b = (long)(rgb - 1000000.0*r - 1000.0*g);

  *mat_color_R = r;
  *mat_color_G = g;
  *mat_color_B = b;

  return 0;
}

/*****************************************************************************/

/*****************************************************************************/

/* Set the x,y,z,time coordinates in the posInfo -struct prior, and the */
/* routine fills  the rest of the fields. Time is ignored. */

long getPositionInfo(long nPositions,pinfo* posInfo)
{
  long err;
  char *uniname;
  char *celname;
  char *matname;

  /*Copy the coordinates over*/

  long iPos;

#pragma omp parallel for private(iPos,uniname,celname,matname) shared(posInfo)

  for(iPos=0; iPos<nPositions; iPos++ )
    {
      err = MatPosUCM(posInfo[iPos].x,
                      posInfo[iPos].y,
                      posInfo[iPos].z,
                      0.0,
                      &(posInfo[iPos].universe),
                      &(posInfo[iPos].cell),
                      &(posInfo[iPos].material),
                      &uniname,&celname,&matname,
                      &(posInfo[iPos].material_color_RGB[0]),
                      &(posInfo[iPos].material_color_RGB[1]),
                      &(posInfo[iPos].material_color_RGB[2])
                      );

      /* The strings need to be copied over to the python data */

      strncpy(posInfo[iPos].universe_name,uniname,MAX_STR);
      posInfo[iPos].universe_name[MAX_STR-1]='\0'; /*Make sure we always terminate the string*/

      strncpy(posInfo[iPos].cell_name,celname,MAX_STR);
      posInfo[iPos].cell_name[MAX_STR-1]='\0'; /*Make sure we always terminate the string*/

      strncpy(posInfo[iPos].material_name,matname,MAX_STR);
      posInfo[iPos].material_name[MAX_STR-1]='\0'; /*Make sure we always terminate the string*/

    }

  return 0;
}

/*****************************************************************************/

/*****************************************************************************/

long getGeometryParameters(struct geom_stats *stats)
{
  long gpl;
  double xmin, xmax, ymin, ymax, zmin, zmax, tmp;

  /* Get pointer to pstructure */

  gpl = (long)RDB[DATA_PTR_GPL0];
  CheckPointer(FUNCTION_NAME, "(gpl)", DATA_ARRAY, gpl);

  /* Get boundaries of geometry */

  xmin = RDB[DATA_GEOM_MINX];
  xmax = RDB[DATA_GEOM_MAXX];
  ymin = RDB[DATA_GEOM_MINY];
  ymax = RDB[DATA_GEOM_MAXY];
  zmin = RDB[DATA_GEOM_MINZ];
  zmax = RDB[DATA_GEOM_MAXZ];

  /* Check if boundaries are set by user */

  if (RDB[gpl + GPL_XMIN] != -INFTY)
    xmin = RDB[gpl + GPL_XMIN];

  if (RDB[gpl + GPL_XMAX] !=  INFTY)
    xmax = RDB[gpl + GPL_XMAX];

  if (RDB[gpl + GPL_YMIN] != -INFTY)
    ymin = RDB[gpl + GPL_YMIN];

  if (RDB[gpl + GPL_YMAX] !=  INFTY)
    ymax = RDB[gpl + GPL_YMAX];

  if (RDB[gpl + GPL_ZMIN] != -INFTY)
    zmin = RDB[gpl + GPL_ZMIN];

  if (RDB[gpl + GPL_ZMAX] !=  INFTY)
    zmax = RDB[gpl + GPL_ZMAX];

  /* Check boundaries and swap */

  if (xmin > xmax)
    {
      tmp = xmax;
      xmax = xmin;
      xmin = tmp;
    }

  if (ymin > ymax)
    {
      tmp = ymax;
      ymax = ymin;
      ymin = tmp;
    }

  if (zmin > zmax)
    {
      tmp = zmax;
      zmax = zmin;
      zmin = tmp;
    }

  /* Move limits away from boundaries to avoid numerical errors in */
  /* some systems. */

  xmin = xmin + EXTRAP_L;
  xmax = xmax - EXTRAP_L;
  ymin = ymin + EXTRAP_L;
  ymax = ymax - EXTRAP_L;
  zmin = zmin + EXTRAP_L;
  zmax = zmax - EXTRAP_L;

  /* Put values into structure */

  stats->xmin = xmin;
  stats->xmax = xmax;
  stats->ymin = ymin;
  stats->ymax = ymax;
  stats->zmin = zmin;
  stats->zmax = zmax;

  /* Exit OK */

  return 0;
}

/******************************************************************************/

long getNumberOfMaterials(long *nMaterials)
{
  long mat;

  /* Pointer to materials */

  mat = (long)RDB[DATA_PTR_M0];
  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

  /* Get list size */

  *nMaterials = ListSize(mat);

  /* Exit OK */

  return 0;
}

/*****************************************************************************/

long getMaterials(struct material *materials )
{
  long nMaterials, status, i, mat, rgb, r, g, b;

  /* Get number of materials */

  if ((status = getNumberOfMaterials(&nMaterials)) != 0)
    return status;

  /* Reset count */

  i = 0;

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Put density */

      materials[i].density = RDB[mat + MATERIAL_ADENS];

      /* Decompose rgb */

      rgb = (long)RDB[mat + MATERIAL_RGB];

      r = (long)(rgb/1000000.0);
      g = (long)((rgb - 1000000.0*r)/1000.0);
      b = (long)(rgb - 1000000.0*r - 1000.0*g);

      /* Put values */

      materials[i].color[0] = r;
      materials[i].color[1] = g;
      materials[i].color[2] = b;

      /* Put name */

      snprintf(materials[i].name, MAX_STR, "%s",
               GetText(mat + MATERIAL_PTR_NAME));

      /* Update index */

      i++;

      /* Pointer to next */

      mat = NextItem(mat);
    }

  /* Exit OK */

  return 0;
}

/*****************************************************************************/

int free_interactive()
{

  FreeMem();
  return 0;
}

/*****************************************************************************/

void init_plot(void)
{
  long ptr;

  /* Check if no geometry plots are defined */

  if ((long)RDB[DATA_PTR_GPL0] < VALID_PTR)
    {
      /* Create a default plot structure */

      ptr = NewItem(DATA_PTR_GPL0, GPL_BLOCK_SIZE);

      /* Init values */

      WDB[ptr + GPL_XMIN] = -INFTY;
      WDB[ptr + GPL_XMAX] =  INFTY;
      WDB[ptr + GPL_YMIN] = -INFTY;
      WDB[ptr + GPL_YMAX] =  INFTY;
      WDB[ptr + GPL_ZMIN] = -INFTY;
      WDB[ptr + GPL_ZMAX] =  INFTY;
      WDB[ptr + GPL_POS] =  -INFTY;
      WDB[ptr + GPL_PLOT_BOUND] = 2.0;
      WDB[ptr + GPL_PIX_X] = 500.0;
      WDB[ptr + GPL_PIX_Y] = 500.0;
      WDB[ptr + GPL_IDX] = 1.0;
      WDB[DATA_N_GEOM_PLOTS] = 1.0;
    }
}

/*****************************************************************************/

int init_interactive(int argc, char** argv)
{
  long ptr, idx[10000], ncoef;
  char str[MAX_STR];
  double t;

  long n;
  unsigned long seed;

  /* Reset coefficient index and transport time */

  ncoef = -1;
  t = 0.0;

  /* Initialise MPI */

  InitMPI(argc, argv);

  /* Init data */

  InitData();

  /* Initialise signal handler */

  InitSignal();

  /***************************************************************************/

  /***** Initial processing before MPI parallelization ***********************/

  /* Check MPI id number */

  if (mpiid == 0)
    {
      /* Get system stat */

      SystemStat();

      /* Process command line */

      if (ParseCommandLine(argc, argv) < 0)
        return 0;

      /* Init OpenMP related stuff */

      InitOMP();

      /* Remove warning message file */

      sprintf(str, "%s.wrn", GetText(DATA_PTR_INPUT_FNAME));
      remove(str);

      /* Read input */

      ReadInput(GetText(DATA_PTR_INPUT_FNAME));
      fprintf(outp, "\n");

      /* Additional initialization */

      init_plot();

      /* Process surfaces */

      ProcessSurfaces();

      /* Reconfigure complement cells */

      ProcessComplementCells();

      /* Remove void cells */

      RemoveVoidCells();

      /* Read STL geometries */

      ReadSTLGeometry();

      /* Check duplicate input definitions */

      CheckDuplicates();

      /* Read pebble bed geometries */

      ReadPBGeometry();

      /* Read unstructured mesh based geometry */

      ReadUMSHGeometry();

      /* Set optimization */

      SetOptimization();

      /* Initialize secondary RNG */

      srand48(parent_seed);

      /* Update memory size */

      WDB[DATA_TOT_MISC_BYTES] = RDB[DATA_TOT_MISC_BYTES] +
        (double)MemCount();

      /* Process burnup material divisors */

      ProcessDivisors();

      /* Processing for MSR calculations */

      ProcessMSR();

      /* Divide burnable zones */

      DivideBurnMat();

      /* Update memory size */

      WDB[DATA_TOT_MAT_BYTES] = RDB[DATA_TOT_MAT_BYTES] +
        (double)MemCount();

      /* Process stochastic geometries */

      ProcessPBGeometry();

      /* Process unstructured mesh based geometries */

      ProcessUMSHGeometry();

      /* Create universes in geometry */

      CreateGeometry();

      /* Count number of zones */

      ZoneCount(-1, -1, 0);

      /* Create super-imposed search meshes */

      ProcessCellMesh();

      /* Process universe transformations */

      ProcessTransformations();

      /* Process universe symmetries */

      ProcessSymmetries();

      /* Check and remove unused stuff */

      CheckUnused();

      /* Process cells */

      ProcessCells();

      /* Process nests */

      ProcessNests();

      /* Process STL geometries */

      ProcessSTLGeometry();

      /* Set material pointers */

      FindMaterialPointers();

      /* Process boundary conditions (must be called after symmetries */
      /* and transformations are processed) */

      ProcessBC();

      /* Update memory size */

      WDB[DATA_TOT_MISC_BYTES] = RDB[DATA_TOT_MISC_BYTES] +
        (double)MemCount();

      /* Make depletion zones */

      MakeDepletionZones(-1, -1, 0, 0, 0, idx, 0.0, 0.0, 0.0, 0);

      /* Update memory size */

      WDB[DATA_TOT_MAT_BYTES] = RDB[DATA_TOT_MAT_BYTES] +
        (double)MemCount();

      /* Process sources, energy grids and detectors (must be done before */
      /* unused cells and materials are removed) */

      ProcessSources();
      ProcessUserEGrids();

      /* Update memory size */

      WDB[DATA_TOT_MISC_BYTES] = RDB[DATA_TOT_MISC_BYTES] +
        (double)MemCount();

      /* Process detectors */

      ProcessDetectors();

      /* Update memory size */

      WDB[DATA_TOT_RES_BYTES] = RDB[DATA_TOT_RES_BYTES] +
        (double)MemCount();

      /* Remove unused materials */

      ptr = (long)RDB[DATA_PTR_M0];
      RemoveFlaggedItems(ptr, MATERIAL_OPTIONS, OPT_USED, NO);

      /* Process lattices */

      ProcessLattices();

      /* Find universe boundaries */

      UniverseBoundaries();

      /* Calculate nest volumes */

      NestVolumes();

      /* Calculate cell volumes */

      CellVolumes();

      /* Count number of cells */

      CellCount(-1, -1, 0, 1);

      /* Calculate material volumes */

      MaterialVolumes();

      /* MGa: Set MPI id's for decomposed materials */

      SetDDIDSimple();

      /* Print geometry data to output file */

      PrintGeometryData();

      /* Update memory size */

      WDB[DATA_TOT_MISC_BYTES] = RDB[DATA_TOT_MISC_BYTES] +
        (double)MemCount();

      /* Process statistics */

      ProcessStats();

      /* Close surface list */

      if ((ptr = (long)RDB[DATA_PTR_S0]) > VALID_PTR)
        CloseList(ptr);

      /* Close universe list */

      if ((ptr = (long)RDB[DATA_PTR_U0]) > VALID_PTR)
        CloseList(ptr);

      /* Close cell list */

      if ((ptr = (long)RDB[DATA_PTR_C0]) > VALID_PTR)
        CloseList(ptr);

      /* Update memory size */

      WDB[DATA_TOT_RES_BYTES] = RDB[DATA_TOT_RES_BYTES] +
        (double)MemCount();

      /* Process multi-physics interfaces */

      ProcessInterface((long)NO);

      /* Prepare adaptive cell search meshes */

      PrepareCellSearchMesh();

      /* Pre-sort cell lists */

      PreSort();

      /* Process variance reduction and response matrix stuff */
      /* (pitää kutsua ennen ExpandPrivateArrays():tä */

      ProcessVR();
      ProcessRMX();

      /* Expand PRIVA, BUF and RES2 arrays for OpenMP parallel */
      /* calculation */

      ExpandPrivateArrays();
    }

  /* Should the following be done in initialization of the code?*/

  /* Loop over OpenMP threads (this is just to avoid errors) */

  for (n = 0; n < (long)RDB[DATA_OMP_MAX_THREADS]; n++)
    {
      /* Init random number sequence */

      seed = ReInitRNG(n);
      SEED[n*RNG_SZ] = seed;
    }

  /* Expand PRIVA, BUF and RES2 arrays for OpenMP parallel calculation */

  ExpandPrivateArrays();

  return 0;
}

/*****************************************************************************/
