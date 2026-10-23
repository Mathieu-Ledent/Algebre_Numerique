#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "adjacences.h"
#include "normales.h"

/* Verifie calculerAdjacences (FLANN ou exhaustive selon la compilation) contre une recherche
   exhaustive : la k-ieme distance doit coincider (les ex aequo peuvent differer), i est exclu. */
static int testerVoisins(int N, int k, const double pts[][2])
{
    Adjacences adj = { N, k, malloc((size_t)N * k * sizeof(int)) };
    calculerAdjacences(pts, &adj);
    double *d = malloc((size_t)N * sizeof(double));
    int ok = 1;
    for (int i = 0; i < N && ok; i++) {
        double dmaxObtenu = 0.0;
        for (int j = 0; j < k; j++) {
            int v = adj.voisins[i * k + j];
            if (v == i || v < 0 || v >= N) { ok = 0; break; }
            for (int jj = 0; jj < j; jj++) if (adj.voisins[i * k + jj] == v) ok = 0;
            double dd = hypot(pts[v][0] - pts[i][0], pts[v][1] - pts[i][1]);
            if (dd > dmaxObtenu) dmaxObtenu = dd;
        }
        /* k-ieme plus petite distance par selection naive. */
        int m = 0;
        for (int j = 0; j < N; j++) if (j != i) d[m++] = hypot(pts[j][0] - pts[i][0], pts[j][1] - pts[i][1]);
        for (int a = 0; a < k; a++)
            for (int b = a + 1; b < m; b++)
                if (d[b] < d[a]) { double t = d[a]; d[a] = d[b]; d[b] = t; }
        if (fabs(dmaxObtenu - d[k - 1]) > 1e-12) ok = 0;
    }
    printf("Voisins N = %d, k = %d : %s\n", N, k, ok ? "OK" : "ECHEC");
    free(adj.voisins); free(d);
    return ok;
}

int main(void)
{
    const double PI = 3.14159265358979323846;
    int tout = 1;
    const int N = 200;
    double (*pts)[2] = malloc((size_t)N * sizeof *pts);
    srand(1);
    for (int i = 0; i < N; i++) {
        pts[i][0] = (double)rand() / RAND_MAX;
        pts[i][1] = (double)rand() / RAND_MAX;
    }
    tout &= testerVoisins(N, 1, (const double (*)[2])pts);
    tout &= testerVoisins(N, 5, (const double (*)[2])pts);
    tout &= testerVoisins(N, 12, (const double (*)[2])pts);

    /* Cercle : les normales doivent etre radiales et, apres orientation, toutes de meme sens. */
    for (int i = 0; i < N; i++) {
        double t = 2 * PI * i / N;
        pts[i][0] = 0.5 + 0.3 * cos(t);
        pts[i][1] = 0.5 + 0.3 * sin(t);
    }
    Adjacences adj = { N, 6, malloc((size_t)N * 6 * sizeof(int)) };
    double (*bary)[2] = malloc((size_t)N * sizeof *bary), (*n)[2] = malloc((size_t)N * sizeof *n);
    calculerAdjacences((const double (*)[2])pts, &adj);
    calculerDroitesMoyennes((const double (*)[2])pts, &adj, bary, n);
    double erreurRadiale = 0.0;
    for (int i = 0; i < N; i++) {
        double rx = (pts[i][0] - 0.5) / 0.3, ry = (pts[i][1] - 0.5) / 0.3;
        double e = 1.0 - fabs(rx * n[i][0] + ry * n[i][1]);
        if (e > erreurRadiale) erreurRadiale = e;
    }
    int ok = erreurRadiale < 1e-3;
    printf("Normales radiales sur un cercle : 1 - |cos| max = %.2e -> %s\n", erreurRadiale, ok ? "OK" : "ECHEC");
    tout &= ok;

    orienterNormales(&adj, n);
    int sortantes = 0;
    for (int i = 0; i < N; i++)
        if ((pts[i][0] - 0.5) * n[i][0] + (pts[i][1] - 0.5) * n[i][1] > 0) sortantes++;
    ok = (sortantes == 0 || sortantes == N);
    printf("Orientation coherente : %d / %d normales sortantes -> %s\n", sortantes, N, ok ? "OK" : "ECHEC");
    tout &= ok;

    free(adj.voisins); free(bary); free(n); free(pts);
    printf("\n%s\n", tout ? "Tous les tests passent." : "Au moins un test echoue.");
    return tout ? 0 : 1;
}
