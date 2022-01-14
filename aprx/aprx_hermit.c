#include "makespl.h"
#include "piv_ge_solver.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

/* UWAGA: liczbę używanych f. bazowych można ustawić przez wartość
          zmiennej środowiskowej APPROX_BASE_SIZE
*/

/*
 * Funkcje bazowe: n - liczba funkcji a,b - granice przedzialu aproksymacji i
 * - numer funkcji x - wspolrzedna dla ktorej obliczana jest wartosc funkcji
 */

double fi(double a, double b, int n, int i, double x){
	double h = (b - a) / (n - 1);
	double h3 = h * h * h;
	int	hi[5] = {i - 2, i - 1, i, i + 1, i + 2};
	double hx[5];

	for (int j = 0; j < 5; j++)
		hx[j] = a + h * hi[j];

	if ((x < hx[0]) || (x > hx[4]))
		return 0;
	else if (x >= hx[0] && x <= hx[1])
		return 1 / h3 * (x - hx[0]) * (x - hx[0]) * (x - hx[0]);
	else if (x > hx[1] && x <= hx[2])
		return 1 / h3 * (h3 + 3 * h * h * (x - hx[1]) + 3 * h * (x - hx[1]) * (x - hx[1]) - 3 * (x - hx[1]) * (x - hx[1]) * (x - hx[1]));
	else if (x > hx[2] && x <= hx[3])
		return 1 / h3 * (h3 + 3 * h * h * (hx[3] - x) + 3 * h * (hx[3] - x) * (hx[3] - x) - 3 * (hx[3] - x) * (hx[3] - x) * (hx[3] - x));
	else			/* if (x > hx[3]) && (x <= hx[4]) */
		return 1 / h3 * (hx[4] - x) * (hx[4] - x) * (hx[4] - x);
}

/* Pierwsza pochodna fi */
double dfi(double a, double b, int n, int i, double x){
	double h = (b - a) / (n - 1);
	double h3 = h * h * h;
	int	hi[5] = {i - 2, i - 1, i, i + 1, i + 2};
	double hx[5];

	for (int j = 0; j < 5; j++)
		hx[j] = a + h * hi[j];

	if ((x < hx[0]) || (x > hx[4]))
		return 0;
	else if (x >= hx[0] && x <= hx[1])
		return 3 / h3 * (x - hx[0]) * (x - hx[0]);
	else if (x > hx[1] && x <= hx[2])
		return 1 / h3 * (3 * h * h + 6 * h * (x - hx[1]) - 9 * (x - hx[1]) * (x - hx[1]));
	else if (x > hx[2] && x <= hx[3])
		return 1 / h3 * (-3 * h * h - 6 * h * (hx[3] - x) + 9 * (hx[3] - x) * (hx[3] - x));
	else			/* if (x > hx[3]) && (x <= hx[4]) */
		return -3 / h3 * (hx[4] - x) * (hx[4] - x);
}

/* Druga pochodna fi */
double d2fi(double a, double b, int n, int i, double x){
	double h = (b - a) / (n - 1);
	double h3 = h * h * h;
	int	hi[5] = {i - 2, i - 1, i, i + 1, i + 2};
	double hx[5];

	for (int j = 0; j < 5; j++)
		hx[j] = a + h * hi[j];

	if ((x < hx[0]) || (x > hx[4]))
		return 0;
	else if (x >= hx[0] && x <= hx[1])
		return 6 / h3 * (x - hx[0]);
	else if (x > hx[1] && x <= hx[2])
		return 1 / h3 * (6 * h - 18 * (x - hx[1]));
	else if (x > hx[2] && x <= hx[3])
		return 1 / h3 * (6 * h  -18 * (hx[3] - x));
	else			/* if (x > hx[3]) && (x <= hx[4]) */
		return 6 / h3 * (hx[4] - x);
}

/* Trzecia pochodna fi */
double d3fi(double a, double b, int n, int i, double x){
	double h = (b - a) / (n - 1);
	double h3 = h * h * h;
	int	hi[5] = {i - 2, i - 1, i, i + 1, i + 2};
	double hx[5];

	for (int j = 0; j < 5; j++)
		hx[j] = a + h * hi[j];

	if ((x < hx[0]) || (x > hx[4]))
		return 0;
	else if (x >= hx[0] && x <= hx[1])
		return 6 / h3;
	else if (x > hx[1] && x <= hx[2])
		return -18 / h3;
	else if (x > hx[2] && x <= hx[3])
		return 18 / h3;
	else			/* if (x > hx[3]) && (x <= hx[4]) */
		return -6 / h3;
}

void make_spl(points_t * pts, spline_t * spl){
	double *x = pts->x;
	double *y = pts->y;
	double a = x[0]; // przypisz wspolrzedna x pierwszego punktu
	double b = x[pts->n - 1]; // przypisz wspolrzedna ostatniego punktu
    int	nb = pts->n - 3 > 10 ? 10 : pts->n - 3;
	char *nbEnv= getenv( "APPROX_BASE_SIZE" );

	if( nbEnv != NULL && atoi( nbEnv ) > 0 )
		nb = atoi( nbEnv );

	matrix_t *eqs = make_matrix(nb, nb + 1); // tworzenie macierzy ukladu rownan

	for (int j = 0; j < nb; j++) {
		for (int i = 0; i < nb; i++) {
			for (int k = 0; k < pts->n; k++) {
				add_to_entry_matrix(eqs, j, i, fi(a, b, nb, i, x[k]) * fi(a, b, nb, j, x[k]));
            }
        }
		for (int k = 0; k < pts->n; k++){
			add_to_entry_matrix(eqs, j, nb, y[k] * fi(a, b, nb, j, x[k]));
        }
	}

	if (piv_ge_solver(eqs)) { // obsluga bledu przy eliminacji gaussa
		spl->n = 0;
		return;
	}

	if (alloc_spl(spl, nb) == 0) {
		for (int i = 0; i < spl->n; i++) {
			double xx = spl->x[i] = a + i*(b-a)/(spl->n-1);
			
			xx+= 10.0 * DBL_EPSILON;  // zabezpieczenie przed ulokowaniem punktu w poprzednim przedziale
			spl->f[i] = 0;
			spl->f1[i] = 0;
			spl->f2[i] = 0;
			spl->f3[i] = 0;
			
			for (int k = 0; k < nb; k++) {
				double ck = get_entry_matrix(eqs, k, nb);
				
				spl->f[i]  += ck * fi  (a, b, nb, k, xx);
				spl->f1[i] += ck * dfi (a, b, nb, k, xx);
				spl->f2[i] += ck * d2fi(a, b, nb, k, xx);
				spl->f3[i] += ck * d3fi(a, b, nb, k, xx);
			}
		}
	}
}
