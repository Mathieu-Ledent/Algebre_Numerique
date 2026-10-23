#ifndef POISSON_H
#define POISSON_H

#include "sparse_matrix.h"
#include "grille.h"




/* Construit la matrice de poisson A, stockée en CSC, avec la numerotation numero[]*/

void construireMatricePoisson(const Grille *grille,
                              const int numero[],
                              MatriceCSC *A);

#endif