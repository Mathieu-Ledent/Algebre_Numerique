#include <stdlib.h>
#include "adjacences.h"

/* Optimisations GCC des noyaux de calcul (ignorees par les autres compilateurs).
   Pas de -ffast-math : on garde l'arithmetique IEEE stricte pour Cholesky. */
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize("O3", "unroll-loops")
#endif

/* Version de reference SANS FLANN : recherche exhaustive en O(N^2) avec selection partielle
   des k plus petites distances au carre (insertion dans une liste triee de taille k).
   Utilisee par le Makefile quand FLANN n'est pas installe, et pour valider la version FLANN. */
void calculerAdjacences(const double points[][2], Adjacences *adj)
{
    const int N = adj->N, k = adj->k;
    double *meilleures = malloc((size_t)k * sizeof(double));
    if (meilleures == NULL) return;

    for (int i = 0; i < N; i++) {
        int *liste = adj->voisins + (size_t)i * k;
        int trouves = 0;
        for (int j = 0; j < N; j++) {
            if (j == i) continue;               /* i est exclu par son indice */
            const double dx = points[j][0] - points[i][0];
            const double dy = points[j][1] - points[i][1];
            const double d2 = dx * dx + dy * dy;
            if (trouves == k && d2 >= meilleures[k - 1]) continue;
            int p = (trouves < k) ? trouves++ : k - 1;
            while (p > 0 && meilleures[p - 1] > d2) {
                meilleures[p] = meilleures[p - 1];
                liste[p] = liste[p - 1];
                p--;
            }
            meilleures[p] = d2;
            liste[p] = j;
        }
    }
    free(meilleures);
}
