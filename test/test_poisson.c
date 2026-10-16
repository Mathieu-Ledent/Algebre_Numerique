#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "poisson.h"

static void afficher(const MatriceCSC *A)
{
    printf("n = %d, nnz = %zu\n", A->n, A->nnz);
    printf("debutColonne = {");
    for (int j = 0; j <= A->n; j++) printf("%zu%s", A->debutColonne[j], j < A->n ? ", " : "}\n");
    printf("lignes       = {");
    for (size_t k = 0; k < A->nnz; k++) printf("%d%s", A->lignes[k], k + 1 < A->nnz ? ", " : "}\n");
    printf("valeurs      = {");
    for (size_t k = 0; k < A->nnz; k++) printf("%g%s", A->valeurs[k], k + 1 < A->nnz ? ", " : "}\n");
}

/* Compare A avec une matrice attendue reconstruite independamment (dense) a partir de la
   geometrie et de numero[]. Verifie aussi les indices de ligne et debutColonne. */
static int verifier(const Grille *g, const int numero[], const MatriceCSC *A)
{
    const int nxi = g->nx - 2, nyi = g->ny - 2, M = nxi * nyi;
    const double invh2 = 1.0 / (g->h * g->h);
    int ok = 1;

    if (A->n != M) { printf("  ECHEC : n = %d au lieu de %d\n", A->n, M); return 0; }

    double *attendu = calloc((size_t)M * M, sizeof(double));
    double *obtenu  = calloc((size_t)M * M, sizeof(double));

    for (int iy = 0; iy < nyi; iy++) {
        for (int ix = 0; ix < nxi; ix++) {
            int p = numero[iy * nxi + ix];
            attendu[(size_t)p * M + p] = 4.0 * invh2;
            if (ix < nxi - 1) {
                int q = numero[iy * nxi + ix + 1];
                int r = p > q ? p : q, c = p > q ? q : p;
                attendu[(size_t)r * M + c] = -invh2;
            }
            if (iy < nyi - 1) {
                int q = numero[(iy + 1) * nxi + ix];
                int r = p > q ? p : q, c = p > q ? q : p;
                attendu[(size_t)r * M + c] = -invh2;
            }
        }
    }

    if (A->debutColonne[0] != 0 || A->debutColonne[M] != A->nnz) {
        printf("  ECHEC : debutColonne[0] ou debutColonne[n] incorrect\n");
        ok = 0;
    }
    for (int c = 0; c < M && ok; c++) {
        for (size_t k = A->debutColonne[c]; k < A->debutColonne[c + 1]; k++) {
            int r = A->lignes[k];
            if (r < c || r >= M) {
                printf("  ECHEC : ligne %d invalide dans la colonne %d\n", r, c);
                ok = 0; break;
            }
            if (k > A->debutColonne[c] && A->lignes[k] <= A->lignes[k - 1]) {
                printf("  ECHEC : lignes non croissantes dans la colonne %d\n", c);
                ok = 0; break;
            }
            obtenu[(size_t)r * M + c] += A->valeurs[k];
        }
        if (ok && A->lignes[A->debutColonne[c]] != c) {
            printf("  ECHEC : la diagonale n'est pas en premier dans la colonne %d\n", c);
            ok = 0;
        }
    }
    for (size_t k = 0; k < (size_t)M * M && ok; k++) {
        if (fabs(attendu[k] - obtenu[k]) > 1e-12 * invh2) {
            printf("  ECHEC : difference a la position (%zu, %zu)\n", k / M, k % M);
            ok = 0;
        }
    }
    free(attendu);
    free(obtenu);
    return ok;
}

/* mode 0 : identite ; 1 : inversee ; 2 : (5 i) % M (permutation si pgcd(5, M) = 1). */
static int tester(const char *nom, int nx, int ny, double h, int mode)
{
    Grille g = { nx, ny, 0.0, 0.0, h, NULL, NULL };
    const int M = (nx - 2) * (ny - 2);
    int *numero = malloc((size_t)M * sizeof(int));
    for (int i = 0; i < M; i++)
        numero[i] = (mode == 0) ? i : (mode == 1) ? M - 1 - i : (5 * i) % M;

    MatriceCSC A;
    construireMatricePoisson(&g, numero, &A);
    int ok = verifier(&g, numero, &A);
    printf("%s : %s\n", nom, ok ? "OK" : "ECHEC");

    libererTableauxCSC(&A);
    free(numero);
    return ok;
}

int main(void)
{
    /* Cas 4x4, h = 1 : a comparer avec la reponse attendue. */
    {
        Grille g = { 4, 4, 0.0, 0.0, 1.0, NULL, NULL };
        int numero[4] = {0, 1, 2, 3};
        MatriceCSC A;
        construireMatricePoisson(&g, numero, &A);
        afficher(&A);
        libererTableauxCSC(&A);
        printf("attendu : nnz = 8, debutColonne = {0, 3, 5, 7, 8}\n");
        printf("          lignes = {0, 1, 2, 1, 3, 2, 3, 3}, valeurs = {4, -1, -1, 4, -1, 4, -1, 4}\n\n");
    }

    int tout = 1;
    tout &= tester("4x4,   identite       ", 4, 4, 1.0, 0);
    tout &= tester("6x5,   h=0.5, identite", 6, 5, 0.5, 0);
    tout &= tester("6x5,   h=0.5, inversee", 6, 5, 0.5, 1);
    tout &= tester("6x5,   h=0.5, (5i)%M  ", 6, 5, 0.5, 2);   /* M = 12, pgcd(5, 12) = 1 */
    tout &= tester("22x22, h=0.05         ", 22, 22, 0.05, 0);
    tout &= tester("3x3 (1 seul inconnu)  ", 3, 3, 1.0, 0);

    printf("\n%s\n", tout ? "Tous les tests passent." : "Au moins un test echoue.");
    return tout ? 0 : 1;
}