/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : writefinixinputfile.c                          */
/*                                                                           */
/* Created:       2015/09/02 (VVa)                                           */
/* Last modified: 2015/09/02 (VVa)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Writes finix_rod.inp, finix_options.inp or                   */
/*              finix_scenario.inp for initialization of a single FINIX rod  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#ifdef FINIX

#define FUNCTION_NAME "WriteFinixInputFile:"

#define FILE_TYPE_ROD      0
#define FILE_TYPE_OPTIONS  1
#define FILE_TYPE_SCENARIO 2

/*****************************************************************************/

void WriteFinixInputFile(long fib, long type)
{
  char infile[MAX_STR], outfile[MAX_STR], name[MAX_STR], line[MAX_STR];
  long flag;
  FILE *fin, *fout;

  /* Create strings for filenames and get block name */
  
  switch (type)
    {
    case FILE_TYPE_ROD:
      sprintf(infile, "%s", GetText(DATA_PTR_FINROD_FNAME));
      sprintf(outfile, "finix_rod.inp");
      sprintf(name, "%s", GetText(fib + FINIX_PTR_RODNAME));

      break;
    case FILE_TYPE_OPTIONS:
      sprintf(infile, "%s", GetText(DATA_PTR_FINOPTI_FNAME));
      sprintf(outfile, "finix_options.inp");
      sprintf(name, "%s", GetText(fib + FINIX_PTR_OPTINAME));

      break;
    case FILE_TYPE_SCENARIO:
      sprintf(infile, "%s", GetText(DATA_PTR_FINSCEN_FNAME));
      sprintf(outfile, "finix_scenario.inp");
      sprintf(name, "%s", GetText(fib + FINIX_PTR_SCENNAME));

      break;
    default:
      Die(FUNCTION_NAME, "Unknown input file type %ld", type);

    }

  /* Open infile for reading */

  if ((fin = fopen(infile,"r")) == NULL)
    Die(FUNCTION_NAME, "Could not open file %s for reading", infile);

  /* Open outfile for writing */

  if ((fout = fopen(outfile,"w")) == NULL)
    Die(FUNCTION_NAME, "Could not open file %s for writing", outfile);

  /* Print "USE <blockname>" definition to the input file */

  fprintf(fout, "USE %s\n", name);

  /* Scan infile for beginning of correct block */

  while (fgets (line, MAX_STR, fin)) 
    {
      if (((strstr(line, "begin ") != NULL) || (strstr(line, "BEGIN ") != NULL)) 
          && (strstr(line, name)))
        {
          fputs(line, fout);
          break;
        }      
    }

  /* Check that we did not reach end of file */

  if (feof(fin))
    Die(FUNCTION_NAME, "Could not find beginning of block \"%s\" from %s", 
        name, infile);

  /* Reset flag */

  flag = 0;

  /* Put block to outfile */

  while (fgets (line, MAX_STR, fin)) 
    {
      /* Put line to outfile */

      fputs(line, fout);

      /* Break if this was the end-of-block line */

      if (((strstr(line, "end ") != NULL) || (strstr(line, "END ") != NULL)) 
          && (strstr(line, name)))
        {
          flag = 1;
          break;
        }
    }

  if (!flag)
    Die(FUNCTION_NAME, "Could not end beginning of block \"%s\" from %s", 
        name, infile);

  /* Close files */

  fclose(fin);
  fclose(fout);

}

#endif
