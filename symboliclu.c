/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : symboliclu.c                                   */
/*                                                                           */
/* Created:       2011/05/03 (MPu)                                           */
/* Last modified: 2015/03/30 (MPu)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Laskee symbolisen LU-hajotelman, jotta Gaussin eliminaatiota */
/*              voidaan soveltaa harvaan ccs-matriisiin                      */
/*              (Fill 2 / Rose & Tarjan 1975)                                */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "SymbolicLU"

/*****************************************************************************/
struct ccsMatrix *SymbolicLU(struct ccsMatrix *cmat){

  /* Palauttaa ccs-matriisin A, joka sisältää myös Gaussissa tulevan fill-inin (E u F) */

  long k,l,i,j,ind,n, m, nnz, rs, len, nnzi; /* len = listan pituus, nnzi = nollasta eroavien lkm rivillä i*/
  long *col, *row, *list, *colA, *rowA, *fill, count; 
  complex *val, *valA; 
  struct ccsMatrix *A; 

  /* Aseta muuttujat matriisien kenttiin */

  n = cmat->n; 
  m = cmat->m; 
  if (n != m)
    Die(FUNCTION_NAME, "burnup matrix not square");
  nnz = cmat->nnz; 

  /* Varaa tila uusille suureille  */

  /* pidä listaa vierekkäisistä solmuista    */
  
  list = (long *)Mem(MEM_ALLOC, n, sizeof(long));

 /* pidä kirjaa jokaisen sarakkeen nz:oista */

  fill = (long *)Mem(MEM_ALLOC, n, sizeof(long));
  
  /* Varaa tilaa joka riviltä (sarakkeesta) fill-iniä varten */
  rs  = 0.1*n;  

/*  rs = 1.0; */
  nnz = nnz + rs*n; /* rs uutta paikkaa joka rivillä   */

  A = ccsMatrixNew(n,n,nnz);
  if (A == NULL)
    Die(FUNCTION_NAME, "Failed to allocate space"); 
  
  row     = cmat->rowind; 
  col     = cmat->colptr;
  val     = cmat->values;

  rowA = A->rowind; /* fill-in tähän matriisiin */
  colA = A->colptr;
  valA = A->values;

  /* Init. */

  for (i=0 ; i < n; i++){ 

    CheckValue(FUNCTION_NAME, "i", "", i, 0, cmat->n - 1);
    CheckValue(FUNCTION_NAME, "i", "", i, 0, A->n);

    nnzi    = col[i+1]-col[i];  /* nz tässä sarakkeessa          */
    colA[i] = col[i] + i*rs;    /* aseta sarakepointterit oikein */

    /* JLE (1.6.2011): mitä jos nnzi == 0? */

    /* (Huomaa, että 1. kierroksella i = 0, eli i*rs = 0... )    */

    /* copy: memcpy(dest, src, size) */
    memcpy(valA + colA[i], val + col[i], nnzi * sizeof(complex));
    memcpy(rowA + colA[i], row + col[i], nnzi * sizeof(long));
    memset(valA + colA[i] + nnzi,  0, rs * sizeof(complex)); 
    memset(rowA + colA[i] + nnzi, -1, rs * sizeof(long));
  }
  colA[n] = nnz; /* aseta vielä viimeisen sarakkeen pointteri */

  /*------------------------------------------------------------------*/
  /* symbolinen hajotelma                                             */
  /*------------------------------------------------------------------*/

  /* Algoritmi perustuu seuraavaan lemmaan: 

     kaari (i,j) kuuluu täytetyn matriisin (A u F) graafiin jos 
     solmusta i on polku solmuun j s.e. kaikki välissä olevan solmut 
     ovat pienempiä kuin min(i,j). 

     Kaikki kaaret (*,j) lasketaan etsimällä ne solmut, joista polku solmuun j s.e. 
     edellinen ehto toteutuu */

  for (i = 1 ; i < n ; i++) { /* aloita toisesta sarakkeesta */

    len = -1; /* tyhjä joukko */

    CheckValue(FUNCTION_NAME, "i", "", i, 0, cmat->n - 1);

    nnzi = col[i+1]-col[i];  /* alkuperäinen nollasta eroavien lkm tässä sarakkeessa */
    
    for (j=col[i] ; j < col[i+1] ; j++){ /* käy läpi nollasta eroavat sarakkeessa i */

    /* 2015/03/30: muutettu A->nnz -1 => A->nnz eli sallii nyt nollasarakkeet matriisin oikeassa reunassa (MPu)*/
    CheckValue(FUNCTION_NAME, "j", "", j, 0, cmat->nnz);

      k       = row[j]; /* a_ki ~= 0 eli k -> i       */ 

      CheckValue(FUNCTION_NAME, "k", "", k, 0, cmat->n - 1);

      fill[k] = 1;      /* otetaan tieto talteen      */
      
      if (k < i){      /* Nyt siis k -> i s.e. k < i  */
        len = len + 1; 

        CheckValue(FUNCTION_NAME, "len", "", len, 0, cmat->n - 1);

        /* lisää k listalle            */

        list[len] = k;
      }
    }

    /*----------------------------------------------------------------*/

    /* Se on tää luuppi johon jumittaa (JLe) */

    count = 0;

    while (len >= 0){   /* Etsitään löytyykö m s.e. m ->... ->k -> i ja k < min(m,i) */

      CheckValue(FUNCTION_NAME, "len", "", len, 0, cmat->n - 1);

      k   = list[len];  /* poista edellinen solmu listalta*/
      len = len - 1;

      /* Testataan ikuinen luuppi (JLe) */

      if (count++ > 100000000)
        Die(FUNCTION_NAME, "Infinite loop?");

      CheckValue(FUNCTION_NAME, "k", "", k, 0, A->n - 1);

      /* Käy läpi listaa vastaavat sarakkeet eli etsi m -> k */
      for (l = colA[k]; l < colA[k+1] ; l++) {

        m = rowA[l]; /* a_mk ~= 0 */

        if (m < 0){
          break;
        }

        if (m > k){
          /* Nyt siis m -> k -> i s.e. k < min(m,i) eli (m,i) lisätään täytettyyn matriisin */
        CheckValue(FUNCTION_NAME, "m", "", m, 0, cmat->n - 1);

          if (fill[m] == 0){ /* Ei ollut siellä valmiiksi */
            /* lisää kaari (m,i)  A:han */

            CheckValue(FUNCTION_NAME, "i", "", i, 0, A->n - 1);

            ind = colA[i] + nnzi;
            if (ind >= colA[i+1]){ /* varaa lisää tilaa */

              /* Tässä nyt uusi funktio */
              ccsMatrixColSpace(A,i,rs); 

              colA = A->colptr;
              rowA = A->rowind;
              valA = A->values;
              nnz  = A->nnz;
            }
            
        /*2015/03/30: Tähän muutettu A->nnz-1 => A->nnz (MPu)*/
            CheckValue(FUNCTION_NAME, "ind", "", ind, 0, A->nnz);

            rowA[ind] = m; /* Nyt matriisissa */
            nnzi = nnzi + 1;

            CheckValue(FUNCTION_NAME, "m", "", m, 0, cmat->n - 1);

            fill[m] = 1;

       /* Seuraava if-osio oli aikaisemmin tämän if-osion ulkopuolella vaikka pitäisi olla sisällä! 
       => selitys päättymättömään luuppiin...? */

        if (m < i){ /* HUOM! edelleen m > k */
            /* Nyt siis m -> k -> i s.e. m & k < i                        */ 
            /* => lisätään m listalle, koska voi löytyä j -> m -> k -> i  */
               len = len + 1;

               CheckValue(FUNCTION_NAME, "len", "", len, 0, cmat->n - 1);

               list[len] = m;
          }
          }
        }                 

        CheckValue(FUNCTION_NAME, "k", "", k, 0, A->n - 1);
      }
    }
    /*-----------------------------------------------------------------*/

    CheckValue(FUNCTION_NAME, "i", "", i, 0, A->n - 1);

    for (l = colA[i] ; l < colA[i] + nnzi ; l++) { /* nollaa fill seuraava kierrosta varten */

      CheckValue(FUNCTION_NAME, "l", "", l, 0, A->nnz - 1);

      k = rowA[l];

      CheckValue(FUNCTION_NAME, "k", "", k, 0, cmat->n - 1);

      fill[k] = 0;

      CheckValue(FUNCTION_NAME, "i", "", i, 0, A->n - 1);
    }
  } /* i silmukka loppuu tähän */
  /*-------------------------------------------------------------------------*/

  ccsMatrixIsort(A);   /* lisäyssorttaus */ 

  Mem(MEM_FREE, list);  
  Mem(MEM_FREE, fill); 
  
  return A;

}
/*****************************************************************************/
