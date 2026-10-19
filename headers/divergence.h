#ifndef GRILLE_H
#define GRILLE_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct {
    int nx, ny; /* Dimensions, bord compris. */
    double xmin, ymin; /* Origine de la grille. */
    double h; /* Pas identique en x et y. */
    double *g; /* Divergence : nx*ny valeurs. */
    double *f; /* Solution : nx*ny valeurs. */
} Grille;

double guassianSplitting(double dist2, double sigma);

void calculerDivergence(int N, const double points[][2],
                        const double normales[][2],
                        const double poids[], double sigma,
                        Grille *grille);

#endif