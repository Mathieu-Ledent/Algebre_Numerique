#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "divergence.h"

static int proche(double a, double b, double tol)
{
    return fabs(a - b) <= tol * fmax(1.0, fmax(fabs(a), fabs(b)));
}

int main(void)
{
    const double pi = 3.14159265358979323846, sigma = 0.4;

    // Validation de la formule de la gaussienne pour des distances connues
    int ok = proche(guassianSplitting(0.0, sigma), 1.0 / (2.0 * pi * sigma * sigma), 1e-14);
    ok &= proche(guassianSplitting(sigma * sigma, sigma),
                 exp(-0.5) / (2.0 * pi * sigma * sigma), 1e-14);

    // Configuration d'un cas test minimal (1 point, 1 normale)
    double point[1][2] = {{0.0, 0.0}};
    double normale[1][2] = {{1.0, 0.0}};
    double poids[1] = {1.0};
    double g1[25], g2[25];
    Grille grille = {5, 5, -1.0, -1.0, 0.5, g1, NULL};

    // Initialisation de la grille avec des NAN pour s'assurer que tout est rempli
    for (int i = 0; i < 25; i++) g1[i] = NAN;
    calculerDivergence(1, point, normale, poids, sigma, &grille);

    // Validation de la structure, des valeurs finies et des propriétés de symétrie/antisymétrie
    for (int s = 0; s < 5; s++) {
        for (int r = 0; r < 5; r++) {
            ok &= isfinite(g1[r + 5 * s]);
            ok &= proche(g1[r + 5 * s], -g1[(4 - r) + 5 * s], 1e-13);
        }
        ok &= proche(g1[2 + 5 * s], 0.0, 1e-13);
    }

    // Validation des poids (poids doublé = divergence doublée)
    for (int i = 0; i < 25; i++) g2[i] = g1[i];
    poids[0] = 2.0;
    calculerDivergence(1, point, normale, poids, sigma, &grille);
    for (int i = 0; i < 25; i++) ok &= proche(g1[i], 2.0 * g2[i], 1e-13);

    // Validation du pointeur de poids NULL (doit retomber sur un poids par défaut de 1.0)
    calculerDivergence(1, point, normale, NULL, sigma, &grille);
    for (int i = 0; i < 25; i++) ok &= proche(g1[i], g2[i], 1e-13);

    // Test de la normalisation discrète et de l'influence du rapport sigma/h
    double erreurGrossiere = 0.0;
    for (int ir = 0; ir < 5; ir++) {
        const double rapport = 0.5 + ir * 0.5;
        const double h = sigma / rapport;
        double somme = 0.0;

        // Sommation de la gaussienne sur un domaine large
        for (double y = -6 * sigma; y <= 6 * sigma + h / 2; y += h)
            for (double x = -6 * sigma; x <= 6 * sigma + h / 2; x += h)
                somme += h * h * guassianSplitting(x * x + y * y, sigma);

        const double erreur = fabs(somme - 1.0);

        // On isole l'erreur pour sigma/h = 0.5 (qui doit être grande)
        if (rapport == 0.5) erreurGrossiere = erreur;
        // Pour des rapports >= 1.0, l'erreur de normalisation doit être très faible
        if (rapport >= 1.0) ok &= erreur < 1e-4;
        printf("sigma/h = %.1f : somme = %.10f, erreur = %.3e\n", rapport, somme, erreur);
    }
    // Vérification que le cas sigma/h = 0.5 présente bien une erreur significative
    ok &= erreurGrossiere > 1e-3;

    printf("Test divergence : %s\n", ok ? "OK" : "ECHEC");
    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
