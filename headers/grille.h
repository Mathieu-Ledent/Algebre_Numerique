#ifndef GRILLE_H
#define GRILLE_H

/* Definition commune de la grille (partagee par divergence, poisson, solver...). */
typedef struct {
    int nx, ny;         /* Dimensions, bord compris. */
    double xmin, ymin;  /* Origine de la grille. */
    double h;           /* Pas identique en x et y. */
    double *g;          /* Divergence : nx*ny valeurs. */
    double *f;          /* Solution : nx*ny valeurs. */
} Grille;

#endif