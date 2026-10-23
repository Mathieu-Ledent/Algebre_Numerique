#ifndef NORMALES_H
#define NORMALES_H
#include "adjacences.h"

void calculerDroitesMoyennes(const double points[][2], const Adjacences *adj, double barycentres[][2], double normales[][2]);

/* Oriente les normales de maniere coherente par parcours en largeur du graphe
   NON oriente des adjacences (i ~ j si j dans N(i) ou i dans N(j)).
   Seul le signe des normales est modifie. */
void orienterNormales(const Adjacences *adj, double normales[][2]);

#endif