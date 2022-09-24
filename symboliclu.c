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

  /* Palauttaa ccs-matriisin A, joka sis�lt�� my�s Gaussissa tulevan fill-inin (E u F) */

  long k,l,i,j,ind,n, m, nnz, rs, len, nnzi; /* len = listan pituus, nnzi = nollasta eroavien lkm rivill� i*/
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

  /* pid� listaa vierekk�isist� solmuista    */
  
  list = (long *)Mem(MEM_ALLOC, n, sizeof(long));

 /* pid� kirjaa jokaisen sarakkeen nz:oista */

  fill = (long *)Mem(MEM_ALLOC, n, sizeof(long));
  
  /* Varaa tilaa joka rivilt� (sarakkeesta) fill-ini� varten */
  rs  = 0.1*n;  

/*  rs = 1.0; */
  nnz = nnz + rs*n; /* rs uutta paikkaa joka rivill�   */

  A = ccsMatrixNew(n,n,nnz);
  if (A == NULL)
    Die(FUNCTION_NAME, "Failed to allocate space"); 
  
  row     = cmat->rowind; 
  col     = cmat->colptr;
  val     = cmat->values;

  rowA = A->rowind; /* fill-in t�h�n matriisiin */
  colA = A->colptr;
  valA = A->values;

  /* Init. */

  for (i=0 ; i < n; i++){ 

    CheckValue(FUNCTION_NAME, "i", "", i, 0, cmat->n - 1);
    CheckValue(FUNCTION_NAME, "i", "", i, 0, A->n);

    nnzi    = col[i+1]-col[i];  /* nz t�ss� sarakkeessa          */
    colA[i] = col[i] + i*rs;    /* aseta sarakepointterit oikein */

    /* JLE (1.6.2011): mit� jos nnzi == 0? */

    /* (Huomaa, ett� 1. kierroksella i = 0, eli i*rs = 0... )    */

    /* copy: memcpy(dest, src, size) */
    memcpy(valA + colA[i], val + col[i], nnzi * sizeof(complex));
    memcpy(rowA + colA[i], row + col[i], nnzi * sizeof(long));
    memset(valA + colA[i] + nnzi,  0, rs * sizeof(complex)); 
    memset(rowA + colA[i] + nnzi, -1, rs * sizeof(long));
  }
  colA[n] = nnz; /* aseta viel� viimeisen sarakkeen pointteri */

  /*------------------------------------------------------------------*/
  /* symbolinen hajotelma                                             */
  /*------------------------------------------------------------------*/

  /* Algoritmi perustuu seuraavaan lemmaan: 

     kaari (i,j) kuuluu t�ytetyn matriisin (A u F) graafiin jos 
     solmusta i on polku solmuun j s.e. kaikki v�liss� olevan solmut 
     ovat pienempi� kuin min(i,j). 

     Kaikki kaaret (*,j) lasketaan etsim�ll� ne solmut, joista polku solmuun j s.e. 
     edellinen ehto toteutuu */

  for (i = 1 ; i < n ; i++) { /* aloita toisesta sarakkeesta */

    len = -1; /* tyhj� joukko */

    CheckValue(FUNCTION_NAME, "i", "", i, 0, cmat->n - 1);

    nnzi = col[i+1]-col[i];  /* alkuper�inen nollasta eroavien lkm t�ss� sarakkeessa */
    
    for (j=col[i] ; j < col[i+1] ; j++){ /* k�y l�pi nollasta eroavat sarakkeessa i */

    /* 2015/03/30: muutettu A->nnz -1 => A->nnz eli sallii nyt nollasarakkeet matriisin oikeassa reunassa (MPu)*/
    CheckValue(FUNCTION_NAME, "j", "", j, 0, cmat->nnz);

      k       = row[j]; /* a_ki ~= 0 eli k -> i       */ 

      CheckValue(FUNCTION_NAME, "k", "", k, 0, cmat->n - 1);

      fill[k] = 1;      /* otetaan tieto talteen      */
      
      if (k < i){      /* Nyt siis k -> i s.e. k < i  */
        len = len + 1; 

        CheckValue(FUNCTION_NAME, "len", "", len, 0, cmat->n - 1);

        /* lis�� k listalle            */

        list[len] = k;
      }
    }

    /*----------------------------------------------------------------*/

    /* Se on t�� luuppi johon jumittaa (JLe) */

    count = 0;

    while (len >= 0){   /* Etsit��n l�ytyyk� m s.e. m ->... ->k -> i ja k < min(m,i) */

      CheckValue(FUNCTION_NAME, "len", "", len, 0, cmat->n - 1);

      k   = list[len];  /* poista edellinen solmu listalta*/
      len = len - 1;

      /* Testataan ikuinen luuppi (JLe) */

      if (count++ > 100000000)
        Die(FUNCTION_NAME, "Infinite loop?");

      CheckValue(FUNCTION_NAME, "k", "", k, 0, A->n - 1);

      /* K�y l�pi listaa vastaavat sarakkeet eli etsi m -> k */
      for (l = colA[k]; l < colA[k+1] ; l++) {

        m = rowA[l]; /* a_mk ~= 0 */

        if (m < 0){
          break;
        }

        if (m > k){
          /* Nyt siis m -> k -> i s.e. k < min(m,i) eli (m,i) lis�t��n t�ytettyyn matriisin */
        CheckValue(FUNCTION_NAME, "m", "", m, 0, cmat->n - 1);

          if (fill[m] == 0){ /* Ei ollut siell� valmiiksi */
            /* lis�� kaari (m,i)  A:han */

            CheckValue(FUNCTION_NAME, "i", "", i, 0, A->n - 1);

            ind = colA[i] + nnzi;
            if (ind >= colA[i+1]){ /* varaa lis�� tilaa */

              /* T�ss� nyt uusi funktio */
              ccsMatrixColSpace(A,i,rs); 

              colA = A->colptr;
              rowA = A->rowind;
              valA = A->values;
              nnz  = A->nnz;
            }
            
        /*2015/03/30: T�h�n muutettu A->nnz-1 => A->nnz (MPu)*/
            CheckValue(FUNCTION_NAME, "ind", "", ind, 0, A->nnz);

            rowA[ind] = m; /* Nyt matriisissa */
            nnzi = nnzi + 1;

            CheckValue(FUNCTION_NAME, "m", "", m, 0, cmat->n - 1);

            fill[m] = 1;

       /* Seuraava if-osio oli aikaisemmin t�m�n if-osion ulkopuolella vaikka pit�isi olla sis�ll�! 
       => selitys p��ttym�tt�m��n luuppiin...? */

        if (m < i){ /* HUOM! edelleen m > k */
            /* Nyt siis m -> k -> i s.e. m & k < i                        */ 
            /* => lis�t��n m listalle, koska voi l�yty� j -> m -> k -> i  */
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
  } /* i silmukka loppuu t�h�n */
  /*-------------------------------------------------------------------------*/

  ccsMatrixIsort(A);   /* lis�yssorttaus */ 

  Mem(MEM_FREE, list);  
  Mem(MEM_FREE, fill); 
  
  return A;

}
/*****************************************************************************/
