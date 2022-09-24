/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : warn.c                                         */
/*                                                                           */
/* Created:       2010/09/14 (JLe)                                           */
/* Last modified: 2017/02/05 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Prints note or warning message for user                      */
/*                                                                           */
/* Comments: - Tähän yhdeksi argumentiksi pointteri laskuriin?               */
/*           - Noi messaget vois kerätä myös erilliseen fileen               */
/*           - Printataan stdouttiin?                                        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "Note:"

/*****************************************************************************/

void Note(long dummy, char * format, ...)
{
  va_list argp;
  va_start (argp, format);

  fprintf(outp, "\n ***** Warning: ");
  vfprintf(outp, format, argp);
  fprintf(outp, "\n\n");
}

/*****************************************************************************/
