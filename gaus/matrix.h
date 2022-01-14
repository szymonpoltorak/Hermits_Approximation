#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <stdio.h>

typedef struct {
	int rn; // wiersze
    int cn; // kolumny
    double *e; // wartosci w macierzy 
} matrix_t;

//tworzy macierz
matrix_t * make_matrix( int rn, int cn );

//czyta macierz
matrix_t * read_matrix( FILE *in );

//wypisz macierz na ekran
void write_matrix( matrix_t *, FILE *out );

//dodaje nowa wartosc do macierzy
void put_entry_matrix( matrix_t *, int i, int j, double val );

//zwieksza wartosc pod danymi indeksami o value
void add_to_entry_matrix( matrix_t *, int i, int j, double val );

double get_entry_matrix( matrix_t *, int i, int j );

//kopiowanie macierzy
matrix_t * copy_matrix( matrix_t *s );

//transponowanie macierzy
matrix_t * transpose_matrix( matrix_t * s );

void xchg_rows( matrix_t *m, int i, int j );

void xchg_cols( matrix_t *m, int i, int j );

matrix_t * mull_matrix( matrix_t *, matrix_t * );

matrix_t * ge_matrix( matrix_t * );

int bs_matrix( matrix_t * );

matrix_t * pivot_ge_matrix( matrix_t *, int *row_per );

void pivot_ge_in_situ_matrix( matrix_t * );

matrix_t * symm_pivot_ge_matrix( matrix_t *, int *per );

int *pivot_get_inv_per( matrix_t *, int *row_per );

#endif
