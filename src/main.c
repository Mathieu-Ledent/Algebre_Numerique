/* Programme principal (section 7.1 de l'enonce), avec quelques options supplementaires :
     -num ligne|rcmk|nd|metis   numerotation des noeuds (defaut : ligne)
     -poids-unitaires           fixe tous les poids a_i a 1 (defaut : poids adaptes)
     -inverser tout|x           inverse toutes les normales ou seulement leur composante x
     -o fichier                 fichier de sortie (defaut : reconstruction.dat)
     -v                         affiche les temps, le remplissage et le residu */
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "grille.h"
#include "reconstruction.h"

int main(int argc, char *argv[])
{
    int k = 5, nGrille = 0, verbeux = 0;
    double sigma = 0.0; /* Zero indique la valeur par defaut. */
    const char *fichier = NULL;
    const char *sortieDat = "reconstruction.dat";
    OptionsReconstruction opt = { NUMEROTATION_LIGNE, 1, 0 };
    int N; /* Nombre de points lu dans le fichier. */
    for (int a = 1; a < argc; ++a) {
        const char *option = argv[a];
        if (strcmp(option, "-N") == 0 ||
            strcmp(option, "-k") == 0 ||
            strcmp(option, "-sigma") == 0) {
            if (++a == argc) goto usage;
            char *finValeur;
            errno = 0;
            if (strcmp(option, "-sigma") == 0) {
                sigma = strtod(argv[a], &finValeur);
                if (errno || finValeur == argv[a] || *finValeur ||
                    !isfinite(sigma) || sigma <= 0.0) goto usage;
            } else {
                long valeur = strtol(argv[a], &finValeur, 10);
                if (errno || finValeur == argv[a] || *finValeur ||
                    valeur < 1 || valeur > INT_MAX) goto usage;
                if (strcmp(option, "-N") == 0)
                    nGrille = (int)valeur;
                else
                    k = (int)valeur;
            }
        } else if (strcmp(option, "-num") == 0) {
            if (++a == argc || lireNumerotation(argv[a], &opt.numerotation) != 0) goto usage;
        } else if (strcmp(option, "-inverser") == 0) {
            if (++a == argc) goto usage;
            if (strcmp(argv[a], "tout") == 0) opt.inverserNormales = 1;
            else if (strcmp(argv[a], "x") == 0) opt.inverserNormales = 2;
            else goto usage;
        } else if (strcmp(option, "-o") == 0) {
            if (++a == argc) goto usage;
            sortieDat = argv[a];
        } else if (strcmp(option, "-poids-unitaires") == 0) {
            opt.poidsAdaptes = 0;
        } else if (strcmp(option, "-v") == 0) {
            verbeux = 1;
        } else {
            if (option[0] == '-' || fichier) goto usage;
            fichier = option;
        }
    }
    if (!fichier || nGrille < 3 ||
        nGrille > INT_MAX / nGrille) goto usage;

    FILE *entree = fopen(fichier, "r");
    if (!entree) {
        perror("Ouverture du nuage");
        return EXIT_FAILURE;
    }
    if (fscanf(entree, "%d", &N) != 1 || N <= k) {
        fprintf(stderr, "Nombre de points invalide.\n");
        fclose(entree);
        return EXIT_FAILURE;
    }
    Grille grille = { .nx = nGrille, .ny = nGrille,
                      .xmin = 0.0, .ymin = 0.0, .h = 0.0,
                      .g = NULL, .f = NULL };
    const size_t taille = (size_t)grille.nx * grille.ny;
    int statut = EXIT_FAILURE;

    /* 1. Allocation des points et du champ de sortie. */
    double (*points)[2] = malloc((size_t)N * sizeof *points);
    grille.f = calloc(taille, sizeof *grille.f);
    if (!points || !grille.f) {
        fprintf(stderr, "Allocation impossible.\n");
        goto fin;
    }

    /* 2. Lecture des N points et de leur boite englobante. */
    double xmin = INFINITY, xmax = -INFINITY;
    double ymin = INFINITY, ymax = -INFINITY;
    for (int i = 0; i < N; ++i) {
        if (fscanf(entree, "%lf %lf",
                   &points[i][0], &points[i][1]) != 2 ||
            !isfinite(points[i][0]) || !isfinite(points[i][1])) {
            fprintf(stderr, "Point %d invalide ou manquant.\n", i);
            goto fin;
        }
        xmin = fmin(xmin, points[i][0]);
        xmax = fmax(xmax, points[i][0]);
        ymin = fmin(ymin, points[i][1]);
        ymax = fmax(ymax, points[i][1]);
    }
    fclose(entree);
    entree = NULL;

    /* Domaine fixe pour comparer plusieurs resolutions. */
    double etendue = fmax(xmax - xmin, ymax - ymin);
    if (!(etendue > 0.0) || !isfinite(etendue)) {
        fprintf(stderr, "Etendue du nuage invalide.\n");
        goto fin;
    }
    grille.h = 1.2 * etendue / (nGrille - 1);
    grille.xmin = xmin - 0.1 * etendue;
    grille.ymin = ymin - 0.1 * etendue;
    if (sigma == 0.0) sigma = 1.5 * grille.h;

    /* 3. Reconstruction complete : champ f et niveau. */
    reglerOptionsReconstruction(&opt);
    double niveau;
    if (reconstructionPoisson(N, (const double (*)[2])points,
                              k, sigma,
                              grille.xmin, grille.ymin, grille.h,
                              grille.nx, grille.ny,
                              grille.f, &niveau) != 0) {
        fprintf(stderr, "Echec de la reconstruction.\n");
        goto fin;
    }

    if (verbeux) {
        const StatistiquesReconstruction *s = statistiquesReconstruction();
        printf("N = %d points, k = %d, sigma = %g, grille %d x %d, h = %g\n",
               N, k, sigma, grille.nx, grille.ny, grille.h);
        printf("Numerotation : %s, poids %s\n", nomNumerotation(opt.numerotation),
               opt.poidsAdaptes ? "adaptes" : "unitaires");
        printf("M = %d inconnues, nnz(A) = %zu, nnz(L) = %zu (remplissage x%.2f)\n",
               s->M, s->nnzA, s->nnzL, (double)s->nnzL / (double)s->nnzA);
        printf("Temps (s) : voisins %.4f | normales %.4f | orientation %.4f | divergence %.4f\n",
               s->tAdjacences, s->tNormales, s->tOrientation, s->tDivergence);
        printf("            numerotation %.4f | assemblage %.4f | symbolique %.4f | numerique %.4f\n",
               s->tNumerotation, s->tAssemblage, s->tSymbolique, s->tNumerique);
        printf("            substitutions %.4f | niveau %.4f | total %.4f\n",
               s->tSubstitutions, s->tNiveau, s->tTotal);
        printf("Residu relatif ||AF - b|| / ||b|| = %.3e\n", s->residuRelatif);
        printf("Niveau tau = %.6e, ecart-type de f_h(p_i) autour de tau = %.3e\n",
               niveau, s->ecartTypeNiveau);
    }

    /* 4. Export du champ et du niveau pour le trace. */
    if (exporterResultats(sortieDat, N,
                          (const double (*)[2])points,
                          &grille, niveau) != 0) {
        fprintf(stderr, "Echec de l'export.\n");
        goto fin;
    }
    statut = EXIT_SUCCESS;

fin:
    /* 5. Liberation, y compris apres une erreur. */
    if (entree) fclose(entree);
    libererCacheReconstruction();
    free(grille.f);
    free(points);
    return statut;

usage:
    fprintf(stderr, "Usage : %s nuage.txt -N taille "
                    "[-sigma longueur] [-k voisins] [-num ligne|rcmk|nd|metis] "
                    "[-poids-unitaires] [-inverser tout|x] [-o sortie.dat] [-v]\n", argv[0]);
    return EXIT_FAILURE;
}
