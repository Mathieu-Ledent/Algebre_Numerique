
#ifndef SOLVER_H
#define SOLVER_H

#include "sparse_matrix.h"
#include "poisson.h"

/* Resout L y = b. y peut etre le meme tableau que b. */
void substitutionAvant(const MatriceCSC *L, const double b[], double y[]);

/* Resout L^T F = y. F peut etre le meme tableau que y. */
void substitutionArriere(const MatriceCSC *L, const double y[], double F[]);

/* Assemble b[numero[i]] = -g au noeud interieur i, resout L L^T F = b,
   puis replace F[numero[i]] dans grille->f (bord mis a 0).
   Renvoie 0 en cas de succes, -1 si les arguments sont invalides ou si une allocation echoue. */
int resoudrePoisson(Grille *grille, const int numero[], const MatriceCSC *L);



#endif