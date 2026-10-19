#ifndef UTILS
#define UTILS

#include "grille.h"
#include "adjacences.h"

/* Temps ecoule (horloge murale) en secondes. */
double chrono(void);

/* Interpolation bilineaire de grille->f au point (x, y).
   Un point hors de la grille est ramene sur son bord. */
double interpolerBilineaire(const Grille *grille, double x, double y);

/* tau = sum a_i f_h(p_i) / sum a_i (poids NULL : a_i = 1). */
double calculerNiveau(const Grille *grille, int N, const double points[][2], const double poids[]);

/* Poids adaptes a l'espacement local : a_i proportionnel a la distance moyenne de p_i a ses
   k voisins, normalise pour que la moyenne des a_i vaille 1. Une zone dense recoit donc
   moins de poids par point. */
void calculerPoidsAdaptes(const double points[][2], const Adjacences *adj, double poids[]);

#endif
