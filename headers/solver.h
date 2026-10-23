
#ifndef SOLVER_H
#define SOLVER_H

#include "sparse_matrix.h"
#include "poisson.h"

/* Resout L y = b. y  */
void substitutionAvant(const MatriceCSC *L, const double b[], double y[]);

/* Resout L^T F = y. */
void substitutionArriere(const MatriceCSC *L, const double y[], double F[]);

/* Utilise le second membre contenu dans grille->g et écrit la solution dans grille->f
Les valeurs de f au bords sont initialisées à 0 et on trouve f grâce à Cholesky creux*/
int resoudrePoisson(Grille *grille, const int numero[], const MatriceCSC *L);



#endif