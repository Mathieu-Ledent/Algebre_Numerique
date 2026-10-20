#include <stdio.h>
#include <stdlib.h>
#include "numerotation.h"
#include "poisson.h"
#include "cholesky.h"

/* numero doit etre une permutation de 0..M-1 ; on affiche le remplissage de L. */
static int tester(int nx, MethodeNumerotation m)
{
    const int nxi = nx - 2, M = nxi * nxi;
    int *numero = malloc((size_t)M * sizeof(int));
    if (calculerNumerotation(m, nxi, nxi, numero) != 0) {
        printf("%-20s nx = %4d : indisponible\n", nomNumerotation(m), nx);
        free(numero);
        return 1;   /* METIS non compile : pas un echec */
    }
    char *vu = calloc((size_t)M, 1);
    int ok = 1;
    for (int i = 0; i < M; i++) {
        if (numero[i] < 0 || numero[i] >= M || vu[numero[i]]) { ok = 0; break; }
        vu[numero[i]] = 1;
    }
    size_t nnzL = 0;
    if (ok) {
        Grille g = { nx, nx, 0.0, 0.0, 1.0, NULL, NULL };
        MatriceCSC A = {0}, L = {0};
        construireMatricePoisson(&g, numero, &A);
        choleskySymbolique(&A, &L);
        ok = (choleskyNumerique(&A, &L) == 0);
        nnzL = L.nnz;
        libererTableauxCSC(&A); libererTableauxCSC(&L);
    }
    printf("%-20s nx = %4d : permutation %s, nnz(L) = %zu\n",
           nomNumerotation(m), nx, ok ? "OK" : "ECHEC", nnzL);
    free(numero); free(vu);
    return ok;
}

int main(void)
{
    int tout = 1;
    const int tailles[] = { 3, 4, 7, 20, 65 };
    for (size_t t = 0; t < sizeof tailles / sizeof *tailles; t++)
        for (int m = NUMEROTATION_LIGNE; m <= NUMEROTATION_METIS; m++)
            tout &= tester(tailles[t], (MethodeNumerotation)m);
    printf("\n%s\n", tout ? "Tous les tests passent." : "Au moins un test echoue.");
    return tout ? 0 : 1;
}
