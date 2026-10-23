#include "sparse_matrix.h"
#include <stddef.h>
#include <stdlib.h>

/* Optimisations GCC des noyaux de calcul (ignorees par les autres compilateurs).
   Pas de -ffast-math : on garde l'arithmetique IEEE stricte pour Cholesky. */
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize("O3", "unroll-loops")
#endif



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

            /* Contribution symétrique, sans recompter l'élément de la diagonale */
            if (i != j) {
                y[j] += a * x[i];
            }
        }
    }
}
