/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : findrowindexes.c                               */
/*                                                                           */
/* Created:       2011/08/17 (JLe)                                           */
/* Last modified: 2012/12/28 (JLe)                                           */
/* Version:       2.1.12                                                     */
/*                                                                           */
/* Description: Finds direct row indexes to speed up NumericGauss()          */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "FindRowIndexes:"

/*****************************************************************************/

void FindRowIndexes(struct ccsMatrix *A)
{
  long n, m, i, j, k, i0, *first;

  /* Varataan muistia */

  A->rowptr = (long *)Mem(MEM_REALLOC, A->rowptr, (A->m + 1)*sizeof(long));
  A->next = (long *)Mem(MEM_REALLOC, A->next, A->nnz*sizeof(long));
  A->colind = (long *)Mem(MEM_REALLOC, A->colind, A->nnz*sizeof(long));
  first = (long *)Mem(MEM_ALLOC, A->n, sizeof(long));

  /* Resetoidaan pointterit ja indeksit */
  
  memset(A->rowptr, 0, (A->m + 1)*sizeof(long));
  memset(A->next, -1, A->nnz*sizeof(long));
  memset(A->colind, 0, A->nnz*sizeof(long));

  /* Haetaan sarakkeiden ensimmäisen nollasta poikkeavan alkion indeksi */

  for (n = 0; n < A->n; n++)
    first[n] = A->colptr[n];

  /* Silmukka rivien yli */

  for (m = 0; m < A->m; m++)
    {
      i0 = -1;

      /* Silmukka sarakkeiden yli */
      
      for (n = 0; n < A->n; n++)
        {
          /* Indeksi sarakedatan alkuun */

          i = first[n];

          /* Tarkistetaan että ollaan vielä samassa sarakkeessa */

          if (((n < A->n - 1) && (i < A->colptr[n + 1])) ||
              ((n == A->n - 1) && (i < A->nnz)))
            {
              /* Ollaanko tällä rivillä? */

              if (A->rowind[i] == m)
                {
                  /* Päivitetään edellisen pointteri seuraavaan ja otetaan */ 
                  /* rivin ensimmäinen pointteri muistiin */

                  if (i0 > -1)
                    A->next[i0]  = i;
                  else
                    A->rowptr[m] = i;

                  /* Sarakkeen indeksi */

                  A->colind[i] = n;

                  /* Päivitetään pointteri edelliseen */
                  
                  i0 = i;
                  
                  /* Päivitetään alkupointteri */
                  
                  first[n]++;
                }
            }
        }
    }

  /* Null pointteri loppuun */

  A->rowptr[m] = -1;

  /* Vapautetaan muisti */

  Mem(MEM_FREE, first);

  /* Nollataan muuttujat ja huijataan kääntäjää olemaan antamatta turhaa */
  /* varoitusta */

  j = 1; 
  j = j - 1;
  k = 1; 
  k = k - 1;

  /* Tää vie älyttömästi aikaa, ei tehdä sitä normaalissa DEBUG-moodissa */

#ifdef DEBUG2

  /* Tarkistetaan */

  for (i = 0; i < A->nnz; i++)
    {
      /* Pointteri seuraavaan */

      if ((j = A->next[i]) > -1)
        {
          /* Tarkistetaan että rivi on sama */
          
          if (A->rowind[i] != A->rowind[j])
            Die(FUNCTION_NAME, "Mismatch in row");

          /* Tarkistetaan että sarake on suurempi */

          if (A->colind[i] >= A->colind[j])
            Die(FUNCTION_NAME, "Mismatch in column");

          /* Tarkistetaan onko välissä alkioita */

          for (k = 0; k < A->nnz; k++)
            if (A->rowind[i] == A->rowind[k])
              if ((A->colind[k] > A->colind[i]) && 
                  (A->colind[k] < A->colind[j]))
                Die(FUNCTION_NAME, "Element missed");
        }
      else
        {
          /* Tarkistetaan että alkio on rivin viimeinen */

          for (k = 0; k < A->nnz; k++)
            if (A->rowind[i] == A->rowind[k])
              if (A->colind[k] > A->colind[i])
                Die(FUNCTION_NAME, "Element missed");
        }
    }

  /* Tarkistetaan että rivien ensimmäiset on ensimmäisiä */

  for (m = 0; m < A->m; m++)
    {
      /* Pointteri ensimmäiseen */

      i = A->rowptr[m];

      /* Tarkistus */

      for (k = 0; k < A->nnz; k++)
        if (A->rowind[i] == A->rowind[k])
          if (A->colind[i] > A->colind[k])
            Die(FUNCTION_NAME, "Mismatch in first element");
    }  

#endif

}

/*****************************************************************************/
