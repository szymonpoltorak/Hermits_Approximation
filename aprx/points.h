#ifndef POINTS_H
#define POINTS_H

#include <stdio.h>

typedef struct {
		int n; //indeksy tablic x i y
		double *x; // tablica wspolrzednych x
		double *y; // tablica wspolrzednych y
} points_t;

// wczytywanie punktow
int read_pts_failed ( FILE* inf, points_t *pts);

void free_points(points_t* pts);

#endif
