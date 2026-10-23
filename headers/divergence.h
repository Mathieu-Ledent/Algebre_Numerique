#ifndef DIVERGENCE_H
#define DIVERGENCE_H

#include "grille.h"

double guassianSplitting(double dist2, double sigma);

void calculerDivergence(int N, const double points[][2],
                        const double normales[][2],
                        const double poids[], double sigma,
                        Grille *grille);

#endif