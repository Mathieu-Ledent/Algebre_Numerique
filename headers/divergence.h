#ifndef DIVERGENCE_H
#define DIVERGENCE_H

#include "grille.h"

/*
Methode permettant de calculer la valeur de G_sigma(d) avec, comme input, la distance au carré et sigma*/
double guassianSplitting(double dist2, double sigma);

/*
Cette methode calcule la divergence en chaque noeud de la grille. 
Elle remplit grille->g qui servira de second membre au problème de Poisson.
*/
void calculerDivergence(int N, const double points[][2],
                        const double normales[][2],
                        const double poids[], double sigma,
                        Grille *grille);

#endif