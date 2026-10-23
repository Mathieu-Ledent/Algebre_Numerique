#include <stdio.h>
#include <stdlib.h>
#include "numerotation.h"
#include "poisson.h"
#include "cholesky.h"
#include "solver.h"
#include "utils.h"

/* Compare les numerotations : remplissage et temps de chaque etape, puis cout d'un second
   membre supplementaire quand la permutation et le facteur sont reutilises. */
int main(int argc, char *argv[])
{
    int nmax = (argc > 1) ? atoi(argv[1]) : 257;
    printf("%6s %-20s %10s %12s %7s %9s %9s %9s %9s %9s\n", "nx", "numerotation", "nnz(A)",
           "nnz(L)", "nnzL/A", "t_perm", "t_assem", "t_symb", "t_num", "t_subst");
    for (int nx = 33; nx <= nmax; nx = 2 * nx - 1) {
        const int nxi = nx - 2, M = nxi * nxi;
        for (int m = NUMEROTATION_LIGNE; m <= NUMEROTATION_METIS; m++) {
            int *numero = malloc((size_t)M * sizeof(int));
            double t0 = chrono();
            if (calculerNumerotation((MethodeNumerotation)m, nxi, nxi, numero) != 0) { free(numero); continue; }
            double t1 = chrono();
            Grille g = { nx, nx, 0.0, 0.0, 1.0 / (nx - 1), calloc((size_t)nx * nx, sizeof(double)),
                         calloc((size_t)nx * nx, sizeof(double)) };
            for (int q = 0; q < nx * nx; q++) g.g[q] = 1.0;
            MatriceCSC A = {0}, L = {0};
            construireMatricePoisson(&g, numero, &A);
            double t2 = chrono();
            choleskySymbolique(&A, &L);
            double t3 = chrono();
            int code = choleskyNumerique(&A, &L);
            double t4 = chrono();
            const int repet = 10;   /* seconds membres supplementaires : seules les substitutions */
            for (int r = 0; r < repet && code == 0; r++) resoudrePoisson(&g, numero, &L);
            double t5 = chrono();
            printf("%6d %-20s %10zu %12zu %7.2f %9.4f %9.4f %9.4f %9.4f %9.5f%s\n", nx,
                   nomNumerotation((MethodeNumerotation)m), A.nnz, L.nnz, (double)L.nnz / A.nnz,
                   t1 - t0, t2 - t1, t3 - t2, t4 - t3, (t5 - t4) / repet, code ? "  ECHEC" : "");
            libererTableauxCSC(&A); libererTableauxCSC(&L);
            free(g.g); free(g.f); free(numero);
        }
    }
    return 0;
}
