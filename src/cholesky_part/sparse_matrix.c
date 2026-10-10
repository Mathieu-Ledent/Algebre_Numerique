
#include "sparse_matrix.h"
#include <stdlib.h>

MatriceCSC *allocate_matrix(int n, size_t nnz)
{
    if (n < 0) {
        return NULL;
    }

    MatriceCSC *mat = malloc(sizeof(MatriceCSC));

    if (mat == NULL) {
        return NULL;
    }

    mat->n = n;
    mat->nnz = nnz;

    mat->valeurs = NULL;
    mat->lignes = NULL;
    mat->debutColonne = NULL;

    if (nnz > 0) {
        mat->valeurs = malloc(nnz * sizeof(double));
        mat->lignes = malloc(nnz * sizeof(int));

        if (mat->valeurs == NULL || mat->lignes == NULL) {
            free(mat->valeurs);
            free(mat->lignes);
            free(mat);
            return NULL;
        }
    }

    mat->debutColonne = calloc((size_t)n + 1,
                               sizeof(size_t));

    if (mat->debutColonne == NULL) {
        free(mat->valeurs);
        free(mat->lignes);
        free(mat);
        return NULL;
    }

    return mat;
}

void free_matrix(MatriceCSC *mat)
{
    if (mat == NULL) {
        return;
    }

    free(mat->valeurs);
    free(mat->lignes);
    free(mat->debutColonne);
    free(mat);
}
