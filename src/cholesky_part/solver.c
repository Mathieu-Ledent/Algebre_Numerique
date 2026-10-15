#include <stdlib.h>
#include "solver.h"

void substitutionAvant(const MatriceCSC *L, const double b[], double y[])
{
    const int M = L->n;
    if (y != b)
        for (int i = 0; i < M; i++) y[i] = b[i];

    for (int i = 0; i < M; i++) {
        const size_t deb = L->debutColonne[i], fin = L->debutColonne[i + 1];
        y[i] /= L->valeurs[deb];                       /* diagonale en premier */
        for (size_t p = deb + 1; p < fin; p++)
            y[L->lignes[p]] -= L->valeurs[p] * y[i];   /* contribution aux lignes suivantes */
    }
}

void substitutionArriere(const MatriceCSC *L, const double y[], double F[])
{
    const int M = L->n;
    for (int i = M - 1; i >= 0; i--) {
        const size_t deb = L->debutColonne[i], fin = L->debutColonne[i + 1];
        double s = y[i];
        for (size_t p = deb + 1; p < fin; p++)
            s -= L->valeurs[p] * F[L->lignes[p]];      /* l_ji * F_j, avec j > i deja calcule */
        F[i] = s / L->valeurs[deb];
    }
}

int resoudrePoisson(Grille *grille, const int numero[], const MatriceCSC *L)
{
    if (grille == NULL || numero == NULL || L == NULL) return -1;
    const int nx = grille->nx, ny = grille->ny;
    const int nxi = nx - 2, nyi = ny - 2, M = nxi * nyi;
    if (nxi <= 0 || nyi <= 0 || L->n != M) return -1;

    double *b = malloc((size_t)M * sizeof(double));
    double *y = malloc((size_t)M * sizeof(double));
    double *F = malloc((size_t)M * sizeof(double));
    if (b == NULL || y == NULL || F == NULL) {
        free(b); free(y); free(F);
        return -1;
    }

    /* Indice du noeud (ix, iy) de la grille complete : iy * nx + ix. */
    for (int iy = 0; iy < nyi; iy++)
        for (int ix = 0; ix < nxi; ix++)
            b[numero[iy * nxi + ix]] = -grille->g[(iy + 1) * nx + (ix + 1)];

    substitutionAvant(L, b, y);
    substitutionArriere(L, y, F);

    for (int k = 0; k < nx * ny; k++) grille->f[k] = 0.0;   /* Dirichlet : f = 0 au bord */
    for (int iy = 0; iy < nyi; iy++)
        for (int ix = 0; ix < nxi; ix++)
            grille->f[(iy + 1) * nx + (ix + 1)] = F[numero[iy * nxi + ix]];

    free(b); free(y); free(F);
    return 0;
}