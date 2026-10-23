#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "reconstruction.h"
#include "adjacences.h"
#include "normales.h"
#include "divergence.h"
#include "poisson.h"
#include "cholesky.h"
#include "solver.h"
#include "utils.h"

static OptionsReconstruction options = { NUMEROTATION_LIGNE, 0, 0 };
static StatistiquesReconstruction stats;

/* Numerotation + matrice + facteur pour la derniere grille : a grille et numerotation fixees,
   seuls le second membre (normales, sigma...) change, on reutilise donc tout. */
static struct {
    int valide;
    int nx, ny;
    double h;
    MethodeNumerotation methode;
    int *numero;
    MatriceCSC A, L;
} cache;

void reglerOptionsReconstruction(const OptionsReconstruction *o) { options = *o; }
const StatistiquesReconstruction *statistiquesReconstruction(void) { return &stats; }

void libererCacheReconstruction(void)
{
    free(cache.numero);
    libererTableauxCSC(&cache.A);
    libererTableauxCSC(&cache.L);
    cache.numero = NULL;
    cache.valide = 0;
}

/* Prepare numero, A et L pour la grille (ou reutilise le cache). Renvoie 0 en cas de succes. */
static int preparerSolveur(const Grille *grille)
{
    if (cache.valide && cache.nx == grille->nx && cache.ny == grille->ny &&
        cache.h == grille->h && cache.methode == options.numerotation) {
        stats.factoriseReutilise = 1;
        return 0;
    }
    libererCacheReconstruction();
    stats.factoriseReutilise = 0;

    const int nxi = grille->nx - 2, nyi = grille->ny - 2, M = nxi * nyi;
    cache.numero = malloc((size_t)M * sizeof(int));
    if (cache.numero == NULL) return -1;

    double t = chrono();
    if (calculerNumerotation(options.numerotation, nxi, nyi, cache.numero) != 0) {
        fprintf(stderr, "Numerotation '%s' indisponible.\n", nomNumerotation(options.numerotation));
        libererCacheReconstruction();
        return -1;
    }
    stats.tNumerotation = chrono() - t;

    t = chrono();
    construireMatricePoisson(grille, cache.numero, &cache.A);
    stats.tAssemblage = chrono() - t;
    if (cache.A.n != M) { libererCacheReconstruction(); return -1; }

    t = chrono();
    choleskySymbolique(&cache.A, &cache.L);
    stats.tSymbolique = chrono() - t;
    if (cache.L.n != M) { libererCacheReconstruction(); return -1; }

    t = chrono();
    int code = choleskyNumerique(&cache.A, &cache.L);
    stats.tNumerique = chrono() - t;
    if (code != 0) {
        fprintf(stderr, "Cholesky : pivot non positif en colonne %d.\n", code - 1);
        libererCacheReconstruction();
        return -1;
    }

    cache.valide = 1;
    cache.nx = grille->nx;
    cache.ny = grille->ny;
    cache.h = grille->h;
    cache.methode = options.numerotation;
    return 0;
}

/* Residu relatif ||A F - b||_2 / ||b||_2 avec la matrice originale (pas le facteur). */
static double residuRelatif(const Grille *grille)
{
    const int nx = grille->nx, nxi = nx - 2, nyi = grille->ny - 2, M = nxi * nyi;
    double *F = malloc((size_t)M * sizeof(double));
    double *b = malloc((size_t)M * sizeof(double));
    double *AF = malloc((size_t)M * sizeof(double));
    double r = NAN;
    if (F && b && AF) {
        for (int iy = 0; iy < nyi; iy++)
            for (int ix = 0; ix < nxi; ix++) {
                const int i = iy * nxi + ix, q = (iy + 1) * nx + (ix + 1);
                F[cache.numero[i]] = grille->f[q];
                b[cache.numero[i]] = -grille->g[q];
            }
        produitSymetrique(&cache.A, F, AF);
        double n1 = 0.0, n2 = 0.0;
        for (int i = 0; i < M; i++) { n1 += (AF[i] - b[i]) * (AF[i] - b[i]); n2 += b[i] * b[i]; }
        r = (n2 > 0.0) ? sqrt(n1 / n2) : sqrt(n1);
    }
    free(F); free(b); free(AF);
    return r;
}

int reconstructionPoisson(int N, const double points[][2],
                          int k, double sigma,
                          double xmin, double ymin, double h,
                          int nx, int ny,
                          double *f, double *niveau)
{
    if (N < 2 || points == NULL || k < 1 || k >= N || !(sigma > 0.0) || !(h > 0.0) ||
        nx < 3 || ny < 3 || f == NULL || niveau == NULL)
        return -1;

    const double tDebut = chrono();
    stats = (StatistiquesReconstruction){0};

    Grille grille = { nx, ny, xmin, ymin, h, NULL, f };
    Adjacences adj = { N, k, NULL };
    double (*barycentres)[2] = malloc((size_t)N * sizeof *barycentres);
    double (*normales)[2] = malloc((size_t)N * sizeof *normales);
    double *poids = malloc((size_t)N * sizeof(double));
    adj.voisins = malloc((size_t)N * k * sizeof(int));
    grille.g = calloc((size_t)nx * ny, sizeof(double));
    int statut = -1;

    if (!barycentres || !normales || !poids || !adj.voisins || !grille.g) {
        fprintf(stderr, "Allocation impossible.\n");
        goto fin;
    }

    /* 1. k plus proches voisins. */
    double t = chrono();
    calculerAdjacences(points, &adj);
    stats.tAdjacences = chrono() - t;

    /* 2. Droites moyennes et normales locales. */
    t = chrono();
    calculerDroitesMoyennes(points, &adj, barycentres, normales);
    stats.tNormales = chrono() - t;

    /* 3. Orientation coherente par parcours en largeur. */
    t = chrono();
    orienterNormales(&adj, normales);
    stats.tOrientation = chrono() - t;
    if (options.inverserNormales == 1)
        for (int i = 0; i < N; i++) { normales[i][0] = -normales[i][0]; normales[i][1] = -normales[i][1]; }
    else if (options.inverserNormales == 2)
        for (int i = 0; i < N; i++) normales[i][0] = -normales[i][0];

    /* 4. Poids puis etalement gaussien : second membre g = div n. */
    for (int i = 0; i < N; i++) poids[i] = 1.0;
    if (options.poidsAdaptes) calculerPoidsAdaptes(points, &adj, poids);
    t = chrono();
    calculerDivergence(N, points, (const double (*)[2])normales, poids, sigma, &grille);
    stats.tDivergence = chrono() - t;

    /* 5. Numerotation, assemblage, Cholesky symbolique et numerique (ou cache). */
    if (preparerSolveur(&grille) != 0) goto fin;
    stats.M = cache.A.n;
    stats.nnzA = cache.A.nnz;
    stats.nnzL = cache.L.nnz;

    /* 6. Substitutions et retour a l'ordre de la grille (bord nul). */
    t = chrono();
    if (resoudrePoisson(&grille, cache.numero, &cache.L) != 0) goto fin;
    stats.tSubstitutions = chrono() - t;
    stats.residuRelatif = residuRelatif(&grille);

    /* 7. Niveau tau : moyenne ponderee de f_h aux points du nuage. */
    t = chrono();
    *niveau = calculerNiveau(&grille, N, points, poids);
    double s2 = 0.0;
    for (int i = 0; i < N; i++) {
        const double e = interpolerBilineaire(&grille, points[i][0], points[i][1]) - *niveau;
        s2 += e * e;
    }
    stats.ecartTypeNiveau = sqrt(s2 / N);
    stats.tNiveau = chrono() - t;
    statut = 0;

fin:
    stats.tTotal = chrono() - tDebut;
    free(barycentres);
    free(normales);
    free(poids);
    free(adj.voisins);
    free(grille.g);
    return statut;
}

int exporterResultats(const char *fichier, int N,
                      const double points[][2],
                      const Grille *grille, double niveau)
{
    FILE *sortie = fopen(fichier, "w");
    if (sortie == NULL) { perror(fichier); return -1; }

    fprintf(sortie, "# nuage\n%d\n", N);
    for (int i = 0; i < N; i++)
        fprintf(sortie, "%.17g %.17g\n", points[i][0], points[i][1]);
    fprintf(sortie, "# grille\n%d %d %.17g %.17g %.17g\n",
            grille->nx, grille->ny, grille->xmin, grille->ymin, grille->h);
    fprintf(sortie, "# niveau\n%.17g\n# f\n", niveau);
    for (int s = 0; s < grille->ny; s++)
        for (int r = 0; r < grille->nx; r++)
            fprintf(sortie, "%.17g%c", grille->f[r + grille->nx * s], r + 1 < grille->nx ? ' ' : '\n');

    return fclose(sortie) == 0 ? 0 : -1;
}
