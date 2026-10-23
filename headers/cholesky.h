#ifndef CHOLESKY_H
#define CHOLESKY_H


#include "sparse_matrix.h"

/* On veut résoudre Af = b, pour ce faire, on utilise la factorisation de Cholesky: A = LL^T.
Cholesky symbolique détermine les positions de L qui doivent être stockées sans calculer les coefficients.
*/
void choleskySymbolique(const MatriceCSC *A, MatriceCSC *L);

/* Cholesky numérique va remplir les positions obtenues avec Cholesky symbolique.
- Renvoie 0 en cas de succes.
- Si le pivot de la colonne j est nul ou negatif : renvoie j + 1 (donc > 0).
- Si les arguments sont invalides ou si une allocation echoue : renvoie -1. */
int choleskyNumerique(const MatriceCSC *A, MatriceCSC *L);


#endif