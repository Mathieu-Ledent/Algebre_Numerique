#include "sparse_matrix.h"
#include <stddef.h>
#include <stdlib.h>



MatriceCSC *allocateMatriceCSC(int n, size_t nnz)
{
    if (n <= 0) return NULL;

    MatriceCSC *m = malloc(sizeof(MatriceCSC));
    if (m == NULL) return NULL;

    m->n = n;
    m->nnz = nnz;

    /* calloc(0, ...) peut renvoyer NULL : on alloue au moins 1 element. */
    size_t taille = (nnz > 0) ? nnz : 1;

    m->valeurs      = calloc(taille, sizeof(double));
    m->lignes       = calloc(taille, sizeof(int));
    m->debutColonne = calloc((size_t)n + 1, sizeof(size_t));

    if (m->valeurs == NULL || m->lignes == NULL || m->debutColonne == NULL) {
        free(m->valeurs);
        free(m->lignes);
        free(m->debutColonne);
        free(m);
        return NULL;
    }

    return m;
}

void freeMatriceCSC(MatriceCSC *m)
{
    if (m == NULL) return;
    free(m->valeurs);
    free(m->lignes);
    free(m->debutColonne);
    free(m);
}