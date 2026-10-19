#ifndef POISSON_H
#define POISSON_H

#include "sparse_matrix.h"
#include "grille.h"




/* Construit A = -Delta_h (triangle inferieur, CSC) pour les noeuds interieurs.
   numero[i] = numero, dans la matrice, du noeud interieur i (ordre ligne par ligne).
   A est fournie vide par l'appelant ; ses tableaux sont alloues ici. */
void construireMatricePoisson(const Grille *grille,
                              const int numero[],
                              MatriceCSC *A);

#endif