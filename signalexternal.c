/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : signalexternal.c                               */
/*                                                                           */
/* Created:       2014/07/18 (VVa)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Signals external program and waits for a return signal       */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SignalExternal:"
#define MAX_TRIES 5
/*****************************************************************************/

void SignalExternal(int sigout)
{
  long ppid;
  int insig, tries;
  FILE *fp;

  /* Signalling only from/to MPI task 0 */

  if(mpiid > 0)
    return;

  if(RDB[DATA_CC_SIG_MODE] == (double)SIG_MODE_POSIX)
    {
      /* If no parent program PID is given, cannot signal */

      if((ppid = (long)RDB[DATA_PPID]) == 0)
        return;

      /* Signal the parent that Serpent is waiting for an updated interface*/

      kill((pid_t)ppid, sigout);

      /* Do not wait for interfaces when sending sigterm */

      if(sigout == SIGTERM)
        return;

      /* Print writing message to terminal */

      fprintf(outp, "Waiting for new interfaces");

      fflush(stdout);

      /* Set waiting flag */

      WDB[DATA_WAITING] = (double)YES;

      /* Wait */
      while((long)RDB[DATA_WAITING]  == YES)
        {
          sleep(1);
          fprintf(outp, ".");
          fflush(stdout);
        }

      fprintf(outp, "\n");
    }
  else if(RDB[DATA_CC_SIG_MODE] == (double)SIG_MODE_FILE)
    {

      /* Open out signal file for writing */

      if ((fp = fopen(GetText(DATA_PTR_COM_OUT), "w")) == NULL)
        Error(-1, "Cannot open out communication file \"%s\" for writing",
              GetText(DATA_PTR_COM_OUT));

      /* Write signal to file */

      fprintf(fp, "%d\n", sigout);

      /* Close out signal file */

      fclose(fp);

      /* Do not wait for interfaces when sending sigterm */

      if(sigout == SIGTERM)
        return;

      /* Print writing message to terminal */

      fprintf(outp, "Waiting for new interfaces");

      fflush(stdout);

      /* Set waiting flag */

      WDB[DATA_WAITING] = (double)YES;

      /* Reset the number of unsuccesful tries */

      tries = 0;

      /* Wait */
      while(RDB[DATA_WAITING] == (double)YES)
        {
          sleep(1);
          fprintf(outp, ".");
          fflush(stdout);

          /* Open infile for reading */

          if ((fp = fopen(GetText(DATA_PTR_COM_IN), "r")) == NULL)
            Error(-1, "Cannot open in communication file \"%s\" for reading",
                  GetText(DATA_PTR_COM_IN));

          if (fscanf(fp, "%d", &insig) == EOF)
            {
              /* Check if we have reached the maximum number of tries */

              if (tries == MAX_TRIES)
                {

                  Die(FUNCTION_NAME,"Could not read signal in from "
                      "communication file %s",
                      GetText(DATA_PTR_COM_IN));
                }
              else
                {

                  /* Close the file */

                  fclose(fp);

                  /* Increase the number of unsuccesful tries */

                  tries++;

                  /* Cycle the loop */

                  continue;
                }
            }
          /* Change waiting flag if SIGUSR1 is received */
          /* = external program has updated interfaces  */

          if (insig == SIGUSR1)
            WDB[DATA_WAITING] = (double)NO;
          else if (insig == SIGUSR2)
            {
              /* Change iteration flag if SIGUSR2 is received               */
              /* = external program has deemed convergence to be sufficient */

              WDB[DATA_WAITING] = (double)NO;

              WDB[DATA_ITERATE] = (double)NO;
            }


          /* Close infile to allow writing by external program */

          fclose(fp);
        }

      /* Signal received, reset signal in infile */

      /* Open infile for writing */

      if ((fp = fopen(GetText(DATA_PTR_COM_IN), "w")) == NULL)
        Error(-1, "Cannot open in communication file \"%s\" for writing",
              GetText(DATA_PTR_COM_IN));

      fprintf(fp, "-1\n");

      fclose(fp);

      fprintf(outp, "\n");
    }

  return;
}
