
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "poisson.h"
#include "cholesky.h"
#include "solver.h"

int main(void)
{
    /* Grille 4 x 4 : seuls les 4 noeuds interieurs sont inconnus. */
    const int nx = 4;
    const int ny = 4;

    double g[16] = {0.0};
    double f[16] = {0.0};

    Grille grille = {
        nx, ny,
        0.0, 0.0,
        1.0,
        g, f
    };

    /* Numerotation des 4 noeuds interieurs, ligne par ligne. */
    int numero[4] = {0, 1, 2, 3};

    /*
     * On choisit la solution exacte :
     * F = [1, 2, 3, 4].
     *
     * Pour A avec h = 1, cela donne b = [-1, 3, 7, 11].
     * Comme le solveur construit b = -g, on initialise g = -b.
     */
    g[5]  =  1.0;
    g[6]  = -3.0;
    g[9]  = -7.0;
    g[10] = -11.0;

    MatriceCSC A = {0};
    MatriceCSC L = {0};

    construireMatricePoisson(&grille, numero, &A);

    if (A.n != 4 || A.nnz != 8 ||
        A.valeurs == NULL || A.lignes == NULL ||
        A.debutColonne == NULL) {
        printf("ECHEC : construction de A\n");
        goto erreur;
    }

    choleskySymbolique(&A, &L);

    if (L.n != 4 || L.valeurs == NULL ||
        L.lignes == NULL || L.debutColonne == NULL) {
        printf("ECHEC : factorisation symbolique\n");
        goto erreur;
    }

    if (choleskyNumerique(&A, &L) != 0) {
        printf("ECHEC : factorisation numerique\n");
        goto erreur;
    }

    if (resoudrePoisson(&grille, numero, &L) != 0) {
        printf("ECHEC : resolution du systeme\n");
        goto erreur;
    }

    /* Solution attendue dans la grille complete, bords inclus. */
    const double attendu[16] = {
        0, 0, 0, 0,
        0, 1, 2, 0,
        0, 3, 4, 0,
        0, 0, 0, 0
    };

    for (int i = 0; i < 16; i++) {
        if (fabs(f[i] - attendu[i]) > 1e-10) {
            printf("ECHEC : f[%d] = %.12g, attendu %.12g\n",
                   i, f[i], attendu[i]);
            goto erreur;
        }
    }

    printf("SUCCES : la solution calculee est correcte.\n");
    printf("Valeurs interieures : %.8f %.8f %.8f %.8f\n",
           f[5], f[6], f[9], f[10]);

    free(A.valeurs);
    free(A.lignes);
    free(A.debutColonne);
    free(L.valeurs);
    free(L.lignes);
    free(L.debutColonne);

    return EXIT_SUCCESS;

erreur:
    free(A.valeurs);
    free(A.lignes);
    free(A.debutColonne);
    free(L.valeurs);
    free(L.lignes);
    free(L.debutColonne);

    return EXIT_FAILURE;
}
