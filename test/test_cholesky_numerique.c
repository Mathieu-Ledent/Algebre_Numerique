#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "poisson.h"
#include "cholesky.h"

/* Matrice dense M x M depuis un CSC triangulaire inferieur ; sym = 1 : on symetrise. */
static double *versDense(const MatriceCSC *C, int sym)
{
    const int M = C->n;
    double *D = calloc((size_t)M * M, sizeof(double));
    for (int c = 0; c < M; c++)
        for (size_t p = C->debutColonne[c]; p < C->debutColonne[c + 1]; p++) {
            int r = C->lignes[p];
            D[(size_t)r * M + c] = C->valeurs[p];
            if (sym) D[(size_t)c * M + r] = C->valeurs[p];
        }
    return D;
}

/* max |A - L L^T| / max |A|. */
static double erreurLLt(const MatriceCSC *A, const MatriceCSC *L)
{
    const int M = A->n;
    double *Ad = versDense(A, 1), *Ld = versDense(L, 0);
    double emax = 0.0, amax = 0.0;
    for (int i = 0; i < M; i++)
        for (int j = 0; j < M; j++) {
            double s = 0.0;
            int kmax = (i < j) ? i : j;
            for (int k = 0; k <= kmax; k++)
                s += Ld[(size_t)i * M + k] * Ld[(size_t)j * M + k];
            double e = fabs(s - Ad[(size_t)i * M + j]);
            if (e > emax) emax = e;
            if (fabs(Ad[(size_t)i * M + j]) > amax) amax = fabs(Ad[(size_t)i * M + j]);
        }
    free(Ad); free(Ld);
    return emax / amax;
}

/* mode 0 : identite ; 1 : inversee ; 2 : (7 i) % M (permutation si pgcd(7, M) = 1). */
static int tester(const char *nom, int nx, int ny, double h, int mode)
{
    Grille g = { nx, ny, 0.0, 0.0, h, NULL, NULL };
    const int M = (nx - 2) * (ny - 2);
    int *numero = malloc((size_t)M * sizeof(int));
    for (int i = 0; i < M; i++)
        numero[i] = (mode == 0) ? i : (mode == 1) ? M - 1 - i : (7 * i) % M;

    MatriceCSC A = {0}, L = {0};
    construireMatricePoisson(&g, numero, &A);
    choleskySymbolique(&A, &L);
    int code = choleskyNumerique(&A, &L);
    double err = (code == 0) ? erreurLLt(&A, &L) : -1.0;
    int ok = (code == 0) && (err < 1e-13);
    printf("%s : M = %d, code = %d, erreur relative |A - LL^T| = %.2e -> %s\n",
           nom, M, code, err, ok ? "OK" : "ECHEC");

    libererTableauxCSC(&A);
    libererTableauxCSC(&L);
    free(numero);
    return ok;
}

int main(void)
{
    int tout = 1;

    /* Cas a la main : A = [[4,-1],[-1,4]] -> L = [[2,0],[-0.5,sqrt(3.75)]]. */
    {
        Grille g = { 4, 3, 0.0, 0.0, 1.0, NULL, NULL };   /* 2 x 1 noeuds interieurs */
        int numero[2] = {0, 1};
        MatriceCSC A = {0}, L = {0};
        construireMatricePoisson(&g, numero, &A);
        choleskySymbolique(&A, &L);
        int code = choleskyNumerique(&A, &L);
        int ok = (code == 0 && L.nnz == 3 &&
                  fabs(L.valeurs[0] - 2.0) < 1e-14 &&
                  fabs(L.valeurs[1] + 0.5) < 1e-14 &&
                  fabs(L.valeurs[2] - sqrt(3.75)) < 1e-14);
        printf("Cas 2x2 a la main : %s\n", ok ? "OK" : "ECHEC");
        tout &= ok;
        libererTableauxCSC(&A);
        libererTableauxCSC(&L);
    }

    /* Matrice non definie positive : [[1,2],[2,1]] -> pivot 1 - 4 = -3 en colonne 1 -> code 2. */
    {
        double v[] = {1, 2, 1};
        int    l[] = {0, 1, 1};
        size_t d[] = {0, 2, 3};
        MatriceCSC A = { 2, 3, v, l, d };
        MatriceCSC L = {0};
        choleskySymbolique(&A, &L);
        int code = choleskyNumerique(&A, &L);
        printf("Matrice non SPD : code = %d (attendu 2) -> %s\n", code, code == 2 ? "OK" : "ECHEC");
        tout &= (code == 2);
        libererTableauxCSC(&L);   /* A pointe vers des tableaux locaux : on ne la libere pas */
    }

    tout &= tester("4x4   ligne par ligne ", 4, 4, 1.0, 0);
    tout &= tester("8x6   ligne par ligne ", 8, 6, 0.5, 0);
    tout &= tester("22x22 ligne par ligne ", 22, 22, 0.05, 0);
    tout &= tester("8x6   inversee        ", 8, 6, 0.5, 1);
    tout &= tester("12x12 permutation 7i  ", 12, 12, 0.1, 2);   /* M = 100, pgcd(7, 100) = 1 */
    tout &= tester("3x3   (1 inconnue)    ", 3, 3, 1.0, 0);

    printf("\n%s\n", tout ? "Tous les tests passent." : "Au moins un test echoue.");
    return tout ? 0 : 1;
}