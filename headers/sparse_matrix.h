
#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <stddef.h>

typedef struct MatriceCSC {
    int n;
    size_t nnz;
    double *valeurs;
    int *lignes;
    size_t *debutColonne;
} MatriceCSC;

/* Alloue la structure ET ses trois tableaux (sur le tas).
   Renvoie NULL en cas d'echec. A liberer avec free_matrix. */
MatriceCSC *allocateMatriceCSC(int n, size_t nnz);

/* Libere les trois tableaux ET la structure elle-meme.
 A n'utiliser que sur une matrice creee par allocateMatriceCSC. */
void freeMatriceCSC(MatriceCSC *m);

/* Libere seulement les trois tableaux (pas la structure) et remet les champs a 0/NULL.
   Pour une structure fournie par l'appelant (variable locale), comme A et L. */
void libererTableauxCSC(MatriceCSC *m);

/* y = A x, avec A symetrique dont seul le triangle inferieur est stocke. */
void produitSymetrique(const MatriceCSC *A, const double x[], double y[]);


#endif


