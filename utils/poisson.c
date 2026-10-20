#include <stdlib.h>
#include "poisson.h"

/* Optimisations GCC des noyaux de calcul (ignorees par les autres compilateurs).
   Pas de -ffast-math : on garde l'arithmetique IEEE stricte pour Cholesky. */
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize("O3", "unroll-loops")
#endif



/* MatriceCSC : voir message precedent. */

/*
 * Noeuds interieurs numerotes ligne par ligne : i = iy * (nx-2) + ix,
 * avec ix dans [0, nx-3] et iy dans [0, ny-3].
 * numero[i] = numero du noeud i dans la matrice (hypothese : permutation).
 * Triangle inferieur stocke par colonnes : colonne c contient les lignes r >= c.
 */
void construireMatricePoisson(const Grille *grille, const int numero[], MatriceCSC *A)
{
    const int nxi = grille->nx - 2;
    const int nyi = grille->ny - 2;
    const double invh2 = 1.0 / (grille->h * grille->h);

    A->n = 0; A->nnz = 0;
    A->valeurs = NULL; A->lignes = NULL; A->debutColonne = NULL;
    if (nxi <= 0 || nyi <= 0) return;

    const int M = nxi * nyi;
    /* Chaque voisinage (paire de noeuds) n'est stocke qu'une fois. */
    const size_t nnz = (size_t)M + (size_t)(nxi - 1) * nyi + (size_t)nxi * (nyi - 1);

    A->valeurs      = malloc(nnz * sizeof(double));
    A->lignes       = malloc(nnz * sizeof(int));
    A->debutColonne = calloc((size_t)M + 1, sizeof(size_t));
    size_t *pos     = calloc((size_t)M, sizeof(size_t)); /* prochaine place libre par colonne */

    if (!A->valeurs || !A->lignes || !A->debutColonne || !pos) {
        free(A->valeurs); free(A->lignes); free(A->debutColonne); free(pos);
        A->valeurs = NULL; A->lignes = NULL; A->debutColonne = NULL;
        return;
    }
    A->n = M;
    A->nnz = nnz;

    /* Passe 1 : taille de chaque colonne (diagonale + voisins de numero superieur). */
    for (int iy = 0; iy < nyi; iy++) {
        for (int ix = 0; ix < nxi; ix++) {
            int p = numero[iy * nxi + ix];
            A->debutColonne[p + 1]++;                       /* diagonale */
            if (ix < nxi - 1) {                             /* voisin de droite */
                int q = numero[iy * nxi + ix + 1];
                A->debutColonne[(p < q ? p : q) + 1]++;
            }
            if (iy < nyi - 1) {                             /* voisin du dessus */
                int q = numero[(iy + 1) * nxi + ix];
                A->debutColonne[(p < q ? p : q) + 1]++;
            }
        }
    }
    for (int c = 0; c < M; c++)
        A->debutColonne[c + 1] += A->debutColonne[c];

    /* La diagonale est toujours en premiere position de sa colonne
       (les autres lignes de la colonne c sont > c). */
    for (int c = 0; c < M; c++) {
        size_t k = A->debutColonne[c];
        A->valeurs[k] = 4.0 * invh2;
        A->lignes[k]  = c;
        pos[c] = k + 1;
    }

    /* Passe 2 : coefficients hors diagonale. */
    for (int iy = 0; iy < nyi; iy++) {
        for (int ix = 0; ix < nxi; ix++) {
            int p = numero[iy * nxi + ix];
            int voisins[2], nv = 0;
            if (ix < nxi - 1) voisins[nv++] = numero[iy * nxi + ix + 1];
            if (iy < nyi - 1) voisins[nv++] = numero[(iy + 1) * nxi + ix];
            for (int v = 0; v < nv; v++) {
                int q = voisins[v];
                int c = (p < q) ? p : q;
                int r = (p < q) ? q : p;
                size_t k = pos[c]++;
                A->valeurs[k] = -invh2;
                A->lignes[k]  = r;
            }
        }
    }

    /* Tri par insertion des lignes dans chaque colonne (au plus 5 elements). */
    for (int c = 0; c < M; c++) {
        size_t deb = A->debutColonne[c], fin = A->debutColonne[c + 1];
        for (size_t a = deb + 1; a < fin; a++) {
            int r = A->lignes[a];
            double v = A->valeurs[a];
            size_t b = a;
            while (b > deb && A->lignes[b - 1] > r) {
                A->lignes[b]   = A->lignes[b - 1];
                A->valeurs[b]  = A->valeurs[b - 1];
                b--;
            }
            A->lignes[b]  = r;
            A->valeurs[b] = v;
        }
    }
    free(pos);
}