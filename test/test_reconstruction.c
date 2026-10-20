#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "reconstruction.h"
#include "utils.h"

/* Chaine complete sur un cercle : le centre et un point exterieur doivent etre de part et
   d'autre de tau, et la ligne de niveau doit passer pres du cercle (rayon mesure sur 360 rayons). */
static int tester(int nPts, int nGrille, MethodeNumerotation m, int poids, double bruit)
{
    const double PI = 3.14159265358979323846, R = 1.0;
    double (*pts)[2] = malloc((size_t)nPts * sizeof *pts);
    srand(3);
    for (int i = 0; i < nPts; i++) {
        double t = 2 * PI * i / nPts;
        double r = R + bruit * (2.0 * rand() / RAND_MAX - 1.0);
        pts[i][0] = r * cos(t);
        pts[i][1] = r * sin(t);
    }
    const double etendue = 2 * R;
    Grille G = { nGrille, nGrille, -R - 0.1 * etendue, -R - 0.1 * etendue,
                 1.2 * etendue / (nGrille - 1), NULL, calloc((size_t)nGrille * nGrille, sizeof(double)) };
    OptionsReconstruction o = { m, poids, 0 };
    reglerOptionsReconstruction(&o);
    double tau;
    int code = reconstructionPoisson(nPts, (const double (*)[2])pts, 8, 1.5 * G.h,
                                     G.xmin, G.ymin, G.h, G.nx, G.ny, G.f, &tau);
    int ok = (code == 0);
    double erreurRayon = 0.0;
    if (ok) {
        double fc = interpolerBilineaire(&G, 0.0, 0.0), fe = interpolerBilineaire(&G, 1.15 * R, 0.0);
        ok = (fc - tau) * (fe - tau) < 0.0;
        /* Sur chaque rayon, premier passage par tau en partant du centre (dichotomie). */
        for (int a = 0; a < 360 && ok; a++) {
            double th = 2 * PI * a / 360, lo = 0.0, hi = 1.15 * R;
            double slo = interpolerBilineaire(&G, 0, 0) - tau;
            for (int it = 0; it < 60; it++) {
                double mid = 0.5 * (lo + hi);
                double s = interpolerBilineaire(&G, mid * cos(th), mid * sin(th)) - tau;
                if ((s > 0) == (slo > 0)) lo = mid; else hi = mid;
            }
            double e = fabs(0.5 * (lo + hi) - R);
            if (e > erreurRayon) erreurRayon = e;
        }
        ok = ok && erreurRayon < 3.0 * G.h + 2 * bruit;
    }
    const StatistiquesReconstruction *s = statistiquesReconstruction();
    printf("Cercle %4d pts, N = %3d, %-20s poids %d bruit %.2f : |r - R| max = %.4f (h = %.4f), "
           "residu %.1e, cache %d -> %s\n", nPts, nGrille, nomNumerotation(m), poids, bruit,
           erreurRayon, G.h, s->residuRelatif, s->factoriseReutilise, ok ? "OK" : "ECHEC");
    free(G.f); free(pts);
    return ok;
}

int main(void)
{
    int tout = 1;
    tout &= tester(200, 41, NUMEROTATION_LIGNE, 0, 0.0);
    tout &= tester(200, 41, NUMEROTATION_LIGNE, 1, 0.0);   /* meme grille : facteur reutilise */
    tout &= tester(400, 81, NUMEROTATION_RCMK, 0, 0.0);
    tout &= tester(400, 81, NUMEROTATION_DISSECTION, 0, 0.0);
    tout &= tester(800, 101, NUMEROTATION_DISSECTION, 0, 0.02);
    libererCacheReconstruction();
    printf("\n%s\n", tout ? "Tous les tests passent." : "Au moins un test echoue.");
    return tout ? 0 : 1;
}
