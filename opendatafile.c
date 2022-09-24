/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : opendatafile.c                                 */
/*                                                                           */
/* Created:       2010/11/22 (JLe)                                           */
/* Last modified: 2011/11/11 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Opens data file for reading                                  */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "OpenDataFile:"

/*****************************************************************************/

FILE *OpenDataFile(long ptr, char *type)
{
  long n;
  FILE *fp;
  char *path, *fname, fname2[MAX_STR];

  /* Check pointer */

  if (ptr < 1)
    Die(FUNCTION_NAME, "Pointer error");

  /* Get file name */

  fname = GetText(ptr);

  /* Try absolute path */

  if ((fp = fopen(fname, "r")) != NULL)
    {
      /* Return file pointer */

      return fp;
    }

  /* Get path from environmental variable */

  if ((path = getenv("SERPENT_DATA")) != NULL)
    {
      /* First try */

      sprintf(fname2, "%s%s", path, fname);

      if ((fp = fopen(fname2, "r")) != NULL)
	{
	  /* Put file name */

	  WDB[ptr] = (double)PutText(fname2);

	  /* Return file pointer */

	  return fp;
	}

      /* Second try */

      sprintf(fname2, "%s/%s", path, fname);

      if ((fp = fopen(fname2, "r")) != NULL)
	{
	  /* Put file name */
	  
	  WDB[ptr] = (double)PutText(fname2);
	  
	  /* Return file pointer */
	  
	  return fp;
	}

      /* Parse file from string */

      for (n = strlen(fname) - 1; n > -1; n--)
	if (fname[n] == '/')
	  break;
      
      /* Check if name was succesfully parsed */

      if ((fname[n] == '/') && (n < (long)strlen(fname) - 2))
	{
	  /* Third try */

	  sprintf(fname2, "%s%s", path, &fname[n + 1]);

	  if ((fp = fopen(fname2, "r")) != NULL)
	    {
	      /* Put file name */
	      
	      WDB[ptr] = (double)PutText(fname2);
	      
	      /* Return file pointer */
	      
	      return fp;
	    }
	 	  
	  /* Fourth try */
	  
	  sprintf(fname2, "%s/%s", path, &fname[n + 1]);
	  
	  if ((fp = fopen(fname2, "r")) != NULL)
	    {
	      /* Put file name */
	      
	      WDB[ptr] = (double)PutText(fname2);
	      
	      /* Return file pointer */
	      
	      return fp;
	    }
	}
    }

  /* Error */

  if (path != NULL)
    Error(0, "Unable to locate %s file \"%s\"\n\nEnvironment variable SERPENT_DATA is set to \"%s\"", 
	  type, fname, path);
  else
    Error(0, "Unable to locate %s file \"%s\"\n\nEnvironment variable SERPENT_DATA not set, file path must be absolute", 
	  type, fname);

  /* Terminate */

  exit(-1);
}

/*****************************************************************************/
