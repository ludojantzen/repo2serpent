/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : printtitle.c                                   */
/*                                                                           */
/* Created:       2010/09/10 (JLe)                                           */
/* Last modified: 2017/10/16 (JLe)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Prints logo and title                                        */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "PrintTitle:"

/*****************************************************************************/

void PrintTitle()
{
  char *path;

  fprintf(outp, "\n");
  fprintf(outp, "  _                   .-=-.           .-=-.          .-==-.       \n");
  fprintf(outp, " { }      __        .' O o '.       .' O o '.       /  -<' )--<   \n");
  fprintf(outp, " { }    .' O'.     / o .-. O \\     / o .-. O \\     /  .---`       \n");
  fprintf(outp, " { }   / .-. o\\   /O  /   \\  o\\   /O  /   \\  o\\   /O /            \n");
  fprintf(outp, "  \\ `-` /   \\ O`-'o  /     \\  O`-'o  /     \\  O`-`o /             \n");
  fprintf(outp, "   `-.-`     '.____.'       `._____.'       `.____.'              \n");
  
  fprintf(outp, "\nSerpent 2 beta\n\n");
  fprintf(outp, "A Continuous-energy Monte Carlo Reactor Physics Burnup ");
  fprintf(outp, "Calculation Code\n\n");

  fprintf(outp, " - Version %s (%s) -- Contact: %s\n\n", CODE_VERSION, CODE_DATE, 
         CODE_AUTHOR);

  fprintf(outp, " - Reference: J. Leppanen, et al. \"The Serpent Monte Carlo code: Status,\n              development and applications in 2013.\" Ann. Nucl. Energy,\n              82 (2015) 142-150.\n\n");

#if defined(__DATE__) && defined(__TIME__)

  fprintf(outp, " - Compiled %s %s\n\n", __DATE__, __TIME__);

#endif

#ifdef MPI
  
  fprintf(outp, " - MPI Parallel calculation mode available\n\n");

#else

  fprintf(outp, " - MPI Parallel calculation mode not available\n\n");

#endif

#ifdef OPEN_MP
  
  fprintf(outp, " - OpenMP Parallel calculation mode available\n\n");

#else

  fprintf(outp, " - OpenMP Parallel calculation mode not available\n\n");

#endif

#ifdef NO_GFX_MODE

  fprintf(outp, " - Geometry and mesh plotting not available\n\n");

#else
  
  fprintf(outp, " - Geometry and mesh plotting available\n\n");

#endif

#ifdef DEBUG

  fprintf(outp, " - Source code compiled in debugger mode\n\n");

#endif 

  if ((path = getenv("SERPENT_DATA")) != NULL)
    fprintf(outp, " - Default data path set to: \"%s\"\n\n", path);
  else
    fprintf(outp, " - Default data path not set\n\n");
}

/*****************************************************************************/
