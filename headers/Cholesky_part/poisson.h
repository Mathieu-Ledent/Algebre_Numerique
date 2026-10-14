#include "sparse_matrix.h"


typedef struct {
int nx, ny; /* Dimensions, bord compris. */
double xmin, ymin; /* Origine de la grille. */
double h; /* Pas identique en x et y. */
double *g; /* Divergence : nx*ny valeurs. */
double *f; /* Solution : nx*ny valeurs. */
} Grille;


void construireMatricePoisson(const Grille *grille,
                              const int numero[],
                              MatriceCSC *A);