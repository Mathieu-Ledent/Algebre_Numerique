#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "poisson.h"
#include "cholesky.h"
#include "solver.h"

#ifndef NMAX
#define NMAX 128   /* plus grande grille testee : n = nx - 1 */
#endif

static const double PI = 3.14159265358979323846;

static double secondes(clock_t a, clock_t b) { return (double)(b - a) / CLOCKS_PER_SEC; }
static int pgcd(int a, int b) { while (b) { int t = a % b; a = b; b = t; } return a; }

/* mode 0 : identite ; 1 : inversee ; 2 : (a*i) % M avec a premier avec M. */
static int *creerNumero(int M, int mode)
{
    int *numero = malloc((size_t)M * sizeof(int));
    int a = 7;
    while (pgcd(a, M) != 1) a++;
    for (int i = 0; i < M; i++)
        numero[i] = (mode == 0) ? i : (mode == 1) ? M - 1 - i : (int)(((long long)a * i) % M);
    return numero;
}

/* Solution fabriquee sur [0,1]^2 avec nx = ny. Renvoie 0 si la resolution a reussi. */
static int essai(int nx, int mode, int afficher, double *eh, double *residu)
{
    const int ny = nx, nxi = nx - 2, M = nxi * (ny - 2);
    const double h = 1.0 / (nx - 1);
    Grille G = { nx, ny, 0.0, 0.0, h,
                 malloc((size_t)nx * ny * sizeof(double)),
                 calloc((size_t)nx * ny, sizeof(double)) };
    int *numero = creerNumero(M, mode);

    for (int iy = 0; iy < ny; iy++)
        for (int ix = 0; ix < nx; ix++)
            G.g[iy * nx + ix] = -2.0 * PI * PI * sin(PI * ix * h) * sin(PI * iy * h);

    MatriceCSC A = {0}, L = {0};
    clock_t t0 = clock();
    construireMatricePoisson(&G, numero, &A);
    clock_t t1 = clock();
    choleskySymbolique(&A, &L);
    clock_t t2 = clock();
    int code = choleskyNumerique(&A, &L);
    clock_t t3 = clock();
    int rc = (code == 0) ? resoudrePoisson(&G, numero, &L) : -1;
    clock_t t4 = clock();

    int echec = (code != 0 || rc != 0);
    if (!echec) {
        /* Erreur de discretisation aux noeuds interieurs. */
        double e = 0.0;
        for (int iy = 1; iy < ny - 1; iy++)
            for (int ix = 1; ix < nx - 1; ix++) {
                double ex = sin(PI * ix * h) * sin(PI * iy * h);
                double d = fabs(G.f[iy * nx + ix] - ex);
                if (d > e) e = d;
            }
        *eh = e;

        /* Residu relatif avec la matrice originale, dans la numerotation de l'assemblage. */
        double *F  = malloc((size_t)M * sizeof(double));
        double *b  = malloc((size_t)M * sizeof(double));
        double *AF = malloc((size_t)M * sizeof(double));
        for (int iy = 0; iy < ny - 2; iy++)
            for (int ix = 0; ix < nxi; ix++) {
                int i = iy * nxi + ix, k = (iy + 1) * nx + (ix + 1);
                F[numero[i]] = G.f[k];
                b[numero[i]] = -G.g[k];
            }
        produitSymetrique(&A, F, AF);
        double n1 = 0.0, n2 = 0.0;
        for (int i = 0; i < M; i++) { n1 += (AF[i] - b[i]) * (AF[i] - b[i]); n2 += b[i] * b[i]; }
        *residu = sqrt(n1 / n2);
        free(F); free(b); free(AF);

        if (afficher) {
            printf("   M = %d, nnz(A) = %zu, nnz(L) = %zu\n", M, A.nnz, L.nnz);
            printf("   temps : assemblage %.4fs | symbolique %.4fs | numerique %.4fs | resolution %.4fs\n",
                   secondes(t0, t1), secondes(t1, t2), secondes(t2, t3), secondes(t3, t4));
        }
    }

    libererTableauxCSC(&A);
    libererTableauxCSC(&L);
    free(G.g); free(G.f); free(numero);
    return echec;
}

/* Test algebrique : b = A F*, on doit retrouver F* sans erreur de discretisation. */
static int testAlgebrique(int nx, int mode)
{
    const int ny = nx, nxi = nx - 2, M = nxi * (ny - 2);
    Grille G = { nx, ny, 0.0, 0.0, 1.0 / (nx - 1), NULL, NULL };
    int *numero = creerNumero(M, mode);
    MatriceCSC A = {0}, L = {0};
    construireMatricePoisson(&G, numero, &A);
    choleskySymbolique(&A, &L);
    int code = choleskyNumerique(&A, &L);

    double *Fs = malloc((size_t)M * sizeof(double)), *b = malloc((size_t)M * sizeof(double));
    double *y  = malloc((size_t)M * sizeof(double)), *F = malloc((size_t)M * sizeof(double));
    for (int i = 0; i < M; i++) Fs[i] = sin(1.0 + 0.37 * i);
    produitSymetrique(&A, Fs, b);
    substitutionAvant(&L, b, y);
    substitutionArriere(&L, y, F);

    double emax = 0.0, fmax = 0.0;
    for (int i = 0; i < M; i++) {
        if (fabs(F[i] - Fs[i]) > emax) emax = fabs(F[i] - Fs[i]);
        if (fabs(Fs[i]) > fmax) fmax = fabs(Fs[i]);
    }
    int ok = (code == 0) && (emax / fmax < 1e-9);
    printf("Test algebrique nx = %d, mode %d : erreur relative = %.2e -> %s\n",
           nx, mode, emax / fmax, ok ? "OK" : "ECHEC");

    free(Fs); free(b); free(y); free(F); free(numero);
    libererTableauxCSC(&A);
    libererTableauxCSC(&L);
    return !ok;
}

int main(void)
{
    int echec = 0;

    printf("=== Solution fabriquee, numerotation ligne par ligne ===\n");
    printf("%6s %10s %12s %12s %8s\n", "n", "h", "e_h", "residu rel.", "ordre");
    double prec = 0.0;
    for (int n = 8; n <= NMAX; n *= 2) {
        double eh = 0.0, res = 0.0;
        if (essai(n + 1, 0, n == NMAX, &eh, &res)) {
            printf("ECHEC pour n = %d\n", n);
            echec = 1;
            continue;
        }
        if (prec > 0.0) printf("%6d %10.5f %12.4e %12.2e %8.2f\n", n, 1.0 / n, eh, res, log2(prec / eh));
        else            printf("%6d %10.5f %12.4e %12.2e %8s\n",  n, 1.0 / n, eh, res, "-");
        prec = eh;
    }

    printf("\n=== La numerotation ne doit pas changer la solution ===\n");
    {
        double e0 = 0, e1 = 0, e2 = 0, r = 0;
        int ech = essai(22, 0, 0, &e0, &r) | essai(22, 1, 0, &e1, &r) | essai(22, 2, 0, &e2, &r);
        int ok = !ech && fabs(e0 - e1) < 1e-12 && fabs(e0 - e2) < 1e-12;
        printf("e_h : %.10e | %.10e | %.10e -> %s\n", e0, e1, e2, ok ? "OK" : "ECHEC");
        echec |= !ok;
    }

    printf("\n");
    echec |= testAlgebrique(12, 0);
    echec |= testAlgebrique(12, 1);
    echec |= testAlgebrique(22, 2);

    printf("\n%s\n", echec ? "Au moins un test echoue." : "Tous les tests passent.");
    return echec;
}