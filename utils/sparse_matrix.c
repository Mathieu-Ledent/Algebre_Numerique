#include "sparse_matrix.h"
#include <stddef.h>
#include <stdlib.h>

/* Optimisations GCC des noyaux de calcul (ignorees par les autres compilateurs).
   Pas de -ffast-math : on garde l'arithmetique IEEE stricte pour Cholesky. */
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize("O3", "unroll-loops")
#endif



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

void libererTableauxCSC(MatriceCSC *m)
{
    if (m == NULL) return;

    free(m->valeurs);
    free(m->lignes);
    free(m->debutColonne);

    m->valeurs = NULL;
    m->lignes = NULL;
    m->debutColonne = NULL;

    m->n = 0;
    m->nnz = 0;
}


void produitSymetrique(const MatriceCSC *A,
                       const double x[],
                       double y[])
{
    if (A == NULL || x == NULL || y == NULL) return;

    /* Initialisation du résultat */
    for (int i = 0; i < A->n; i++) {
        y[i] = 0.0;
    }

    /* Parcours des colonnes de la matrice */
    for (int j = 0; j < A->n; j++) {
        for (size_t k = A->debutColonne[j];
             k < A->debutColonne[j + 1];
             k++) {

            int i = A->lignes[k];
            double a = A->valeurs[k];

            /* Contribution du coefficient stocké */
            y[i] += a * x[j];

            /* Contribution symétrique, sauf sur la diagonale */
            if (i != j) {
                y[j] += a * x[i];
            }
        }
    }
}