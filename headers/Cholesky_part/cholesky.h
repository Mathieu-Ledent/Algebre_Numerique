
#include "sparse_matrix.h"

void choleskySymbolique(const MatriceCSC *A, MatriceCSC *L);

/* Renvoie 0 en cas de succes.
   Si le pivot de la colonne j est nul ou negatif : renvoie j + 1 (donc > 0).
   Si les arguments sont invalides ou si une allocation echoue : renvoie -1. */
int choleskyNumerique(const MatriceCSC *A, MatriceCSC *L);