/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : systemstat.c                                   */
/*                                                                           */
/* Created:       2011/05/06 (JLe)                                           */
/* Last modified: 2015/05/16 (JLe)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Reads some system data                                       */
/*                                                                           */
/* Comments: - From Serpent 1.1.9, the values are not really used for        */
/*             anything                                                      */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SystemStat:"

/*****************************************************************************/

void SystemStat()
{
  double sz;
  long n;
  FILE *fp;
  char tmpstr[MAX_STR], unit[MAX_STR];

  /* Try to open /proc/cpuinfo for reading */

  if ((fp = fopen("/proc/cpuinfo", "r")) != NULL)
    {
      /* Skip first 4 lines */
      
      if (fgets(tmpstr, 82, fp) == NULL)
        Die (FUNCTION_NAME, "fgets error");

      if (fgets(tmpstr, 82, fp) == NULL)
        Die (FUNCTION_NAME, "fgets error");

      if (fgets(tmpstr, 82, fp) == NULL)
        Die (FUNCTION_NAME, "fgets error");
      
      if (fgets(tmpstr, 82, fp) == NULL)
        Die (FUNCTION_NAME, "fgets error");

      /* Get model name */

      if (fgets(tmpstr, 82, fp) == NULL)
        Die (FUNCTION_NAME, "fgets error");

      /* Find last character */

      for (n = strlen(tmpstr); n > 0; n--)
        if (isalnum(tmpstr[n]))
          {
            tmpstr[n + 1] = '\0';

            break;
          }

      /* Find first character */

      for (n = 0; n < (int)strlen(tmpstr); n++)
        if (tmpstr[n] == ':')
          break;

      /* Store value */
      
      WDB[DATA_PTR_CPU_NAME] = (double)PutText(&tmpstr[n + 2]);

      /* Skip line */

      if (fgets(tmpstr, 82, fp) == NULL)
        Die (FUNCTION_NAME, "fgets error");

      /* Get MHZ */

      if (fgets(tmpstr, 82, fp) == NULL)
        Die (FUNCTION_NAME, "fgets error");

      /* Find first character */

      for (n = 0; n < (int)strlen(tmpstr); n++)
        if (tmpstr[n] == ':')
          break;

      /* Store value */
      
      WDB[DATA_CPU_MHZ] = atof(&tmpstr[n + 2]);
      
      /* Close file */

      fclose(fp);
    }
  else
    WDB[DATA_PTR_CPU_NAME] = (double)PutText("Unknown");

  /* Try to open /proc/meminfo for reading */

  if ((fp = fopen("/proc/meminfo", "r")) != NULL)
    {
      /* Get total memory */

      if (fscanf(fp, "%s %lf %s\n", tmpstr, &sz, unit) == EOF)
        Die(FUNCTION_NAME, "fscanf error");
      
      /* Check unit */

      if (!strcmp(unit, "kB"))
        WDB[DATA_CPU_MEM] = sz*KILO/GIGA;

      /* Close file */

      fclose(fp);
    }

  /* Working directory */
  
  if (getcwd(tmpstr, sizeof(tmpstr)) == NULL)
    sprintf(tmpstr, "N/A");

  WDB[DATA_PTR_WORKDIR] = (double)PutText(tmpstr);

  /* Host name */
  
  if (gethostname(tmpstr, sizeof(tmpstr)) == -1)
    sprintf(tmpstr, "N/A");

  WDB[DATA_PTR_HOSTNAME] = (double)PutText(tmpstr);

  /* Compile date */

#if defined(__DATE__) && defined(__TIME__)

  sprintf(tmpstr, "%s %s", __DATE__, __TIME__);

#else

  sprintf(tmpstr, "N/A");

#endif

  WDB[DATA_PTR_COMPILE_DATE] = (double)PutText(tmpstr);
}

/*****************************************************************************/
