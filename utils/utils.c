#define _POSIX_C_SOURCE 199309L
#include <math.h>
#include <time.h>
#include "utils.h"

double chrono(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec + 1e-9 * (double)t.tv_nsec;
}

double interpolerBilineaire(const Grille *grille, double x, double y)
{
    const int nx = grille->nx, ny = grille->ny;
    double u = (x - grille->xmin) / grille->h;
    double v = (y - grille->ymin) / grille->h;
    if (u < 0.0) u = 0.0;
    if (v < 0.0) v = 0.0;
    if (u > nx - 1) u = nx - 1;
    if (v > ny - 1) v = ny - 1;

    int r = (int)floor(u), s = (int)floor(v);
    if (r > nx - 2) r = nx - 2;
    if (s > ny - 2) s = ny - 2;
    const double a = u - r, b = v - s;
    const double *f = grille->f;
    return (1 - a) * (1 - b) * f[r + nx * s] + a * (1 - b) * f[r + 1 + nx * s]
         + (1 - a) * b * f[r + nx * (s + 1)] + a * b * f[r + 1 + nx * (s + 1)];
}

double calculerNiveau(const Grille *grille, int N, const double points[][2], const double poids[])
{
    double num = 0.0, den = 0.0;
    for (int i = 0; i < N; i++) {
        const double a = (poids != NULL) ? poids[i] : 1.0;
        num += a * interpolerBilineaire(grille, points[i][0], points[i][1]);
        den += a;
    }
    return num / den;
}

void calculerPoidsAdaptes(const double points[][2], const Adjacences *adj, double poids[])
{
    const int N = adj->N, k = adj->k;
    double somme = 0.0;
    for (int i = 0; i < N; i++) {
        double d = 0.0;
        for (int j = 0; j < k; j++) {
            const int v = adj->voisins[i * k + j];
            d += hypot(points[v][0] - points[i][0], points[v][1] - points[i][1]);
        }
        poids[i] = d / k;
        somme += poids[i];
    }
    const double moyenne = somme / N;
    for (int i = 0; i < N; i++)
        poids[i] = (moyenne > 0.0 && poids[i] > 0.0) ? poids[i] / moyenne : 1.0;
}
