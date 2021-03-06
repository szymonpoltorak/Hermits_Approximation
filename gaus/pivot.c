#include "matrix.h"
#include <stdlib.h>
#include <math.h>

matrix_t *pivot_ge_matrix (matrix_t * a, int *row_per){
  matrix_t *c = copy_matrix (a);
  if (c != NULL) {
    for (int i = 0; i < c->rn; i++)
      row_per[i] = i;
    for (int k = 0; k < c->rn - 1; k++) {      /* eliminujemy (zerujemy) kolumnę nr k */
      int piv = k;              /* wybór eleemntu dominującego - maks. z k-tej kol., poniżej diag */
      
      for (int i = k + 1; i < c->rn; i++)
        if (fabs (*(c->e + i * c->cn + k)) > fabs (*(c->e + piv * c->cn + k)))
          piv = i;
      
      if (piv != k) {           /* jeśli diag. nie jest pivtem - wymień wiersze */
        int tmp;
        xchg_rows (c, piv, k);
        tmp = row_per[k];
        row_per[k] = row_per[piv];
        row_per[piv] = tmp;
      }
      
      for (int i = k + 1; i < c->rn; i++) {    /* pętla po kolejnych
                                           wierszach poniżej diagonalii k,k */
        double d = *(c->e + i * c->cn + k) / *(c->e + k * c->cn + k);
        for (int j = k; j < c->cn; j++)
          *(c->e + i * c->cn + j) -= d * *(c->e + k * c->cn + j);
      }
    }
  }
  return c;
}

void pivot_ge_in_situ_matrix (matrix_t * c){
  for (int k = 0; k < c->rn - 1; k++) {        /* eliminujemy (zerujemy) kolumnę nr k */
    int piv = k;                /* wybór elementu dominującego - maks. z k-tej kol., poniżej diag */
    
    for (int i = k + 1; i < c->rn; i++)
      if (fabs (*(c->e + i * c->cn + k)) > fabs (*(c->e + piv * c->cn + k)))
        piv = i;
    
    if (piv != k) {             /* jeśli diag. nie jest pivtem - wymień wiersze */
      xchg_rows (c, piv, k);
    }
    
    for (int i = k + 1; i < c->rn; i++) {      /* pętla po kolejnych
                                           wierszach poniżej diagonalii k,k */
      double d = *(c->e + i * c->cn + k) / *(c->e + k * c->cn + k);

      for (int j = k; j < c->cn; j++)
        *(c->e + i * c->cn + j) -= d * *(c->e + k * c->cn + j);
    }
  }
}

matrix_t *symm_pivot_ge_matrix (matrix_t * a, int *row_per){
  matrix_t *c = copy_matrix (a);

  if (c != NULL) {
    for (int i = 0; i < c->rn; i++)
      row_per[i] = i;
    
    for (int k = 0; k < c->rn - 1; k++) {      /* eliminujemy (zerujemy) kolumnę nr k */
      int piv = k;              /* wybór eleemntu dominującego - maks. z k-tej kol., poniżej diag */
      
      for (int i = k + 1; i < c->rn; i++)
        if (fabs (*(c->e + i * c->cn + k)) > fabs (*(c->e + piv * c->cn + k)))
          piv = i;
      
      if (piv != k) {           /* jeśli diag. nie jest pivtem - wymień wiersze */
        int tmp;
        xchg_rows (c, piv, k);
        xchg_cols (c, piv, k);
        tmp = row_per[k];
        row_per[k] = row_per[piv];
        row_per[piv] = tmp;
      }

      for (int i = k + 1; i < c->rn; i++) {    /* pętla po kolejnych
                                           wierszach poniżej diagonalii k,k */
        double d = *(c->e + i * c->cn + k) / *(c->e + k * c->cn + k);
        for (int j = k; j < c->cn; j++)
          *(c->e + i * c->cn+ j) -= d * *(c->e + k * c->cn + j);
      }
    }
  }
  return c;
}

int *pivot_get_inv_per (matrix_t * m, int *row_per){
  /* odtwarzamy oryginalną kolejność niewiadowmych */
  int *iper = (int*) malloc (m->rn * sizeof *iper);

  for (int i = 0; i < m->rn; i++)
    iper[row_per[i]] = i;

  return iper;
}
