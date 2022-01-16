#include "makespl.h"
#include "piv_ge_solver.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

double Hn(double x, int n){			
	if ( n == 0 ) {
		return 	1;
	} 
	else if ( n == 1 ) {
		return 	2 * x;
	} 
	else {
        return  2 * x * Hn(x, n - 1) - 2 * Hn(x, n - 1) * Hn(x, n - 2);
	} 
}

double dHn(double x, int n){
	if (n == 0){
		return 0;
	}
	else if (n == 1){
		return 2;
	}
	else {
		return 2 * Hn(x, n - 1) + 2 * x * dHn(x, n - 1) - 2 * (n - 1) * dHn(x, n - 2);
	}
}

double d2Hn(double x, int n){
	if (n == 0){
		return 0;
	}
	else if (n == 1){
		return 0;
	}
	else{
		return 4 * dHn(x, n - 1) + 2 * x * d2Hn(x, n - 1) - 2 * (n - 1) * d2Hn(x, n - 2);
	}
}

double d3Hn(double x, int n){
	if (n == 0){
		return 0;
	}
	else if (n == 1){
		return 0;
	}
	else {
		return 6 * d2Hn(x, n - 1) + 2 * x * d3Hn(x, n - 1) - 2 * (n - 1) * d3Hn(x, n - 2);
	}
}

void make_spl(points_t * pts, spline_t * spl){
	double *x = pts->x;
	double *y = pts->y;
	double a = x[0]; // przypisz wspolrzedna x pierwszego punktu
	double b = x[pts->n - 1]; // przypisz wspolrzedna ostatniego punktu
    int	nb = pts->n - 3 > 10 ? 10 : pts->n - 3;
	char *nbEnv= getenv( "APPROX_BASE_SIZE" );

	if( nbEnv != NULL && atoi( nbEnv ) > 0 ) {
		nb = atoi( nbEnv );
    }

	matrix_t *eqs = make_matrix(nb, nb + 1); // tworzenie macierzy ukladu rownan

	for (int j = 0; j < nb; j++) {
		for (int i = 0; i < nb; i++) {
			for (int k = 0; k < pts->n; k++) {
				add_to_entry_matrix(eqs, j, i, Hn(x[k],i) * Hn(x[k],j));
            }
        }
		for (int k = 0; k < pts->n; k++){
			add_to_entry_matrix(eqs, j, nb, y[k] * Hn(x[k], j));
        }
	}

	if (piv_ge_solver(eqs)) { // obsluga bledu przy eliminacji gaussa
		spl->n = 0;
		return;
	}

	if (alloc_spl(spl, nb) == 0) {
		for (int i = 0; i < spl->n; i++) {
			double xx = spl->x[i] = a + i*(b-a)/(spl->n-1);
			
			xx+= 10.0 * DBL_EPSILON;
			spl->f[i] = 0;
			spl->f1[i] = 0;
			spl->f2[i] = 0;
			spl->f3[i] = 0;
			
			for (int k = 0; k < nb; k++) {
				double ck = get_entry_matrix(eqs, k, nb);
				
				spl->f[i]  += ck * Hn  (xx, k);
				spl->f1[i] += ck * dHn (xx, k);
				spl->f2[i] += ck * d2Hn(xx, k);
				spl->f3[i] += ck * d3Hn(xx, k);
			}
		}
	}
    free(eqs -> e);
    free(eqs);
}
