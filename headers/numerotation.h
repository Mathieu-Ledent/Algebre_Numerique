#ifndef NUMEROTATION_H
#define NUMEROTATION_H

#include "sparse_matrix.h"

/* Methodes de numerotation des noeuds interieurs.
   Dans tous les cas, numero[i] est le numero dans la matrice du noeud interieur i
   (i = iy * nxi + ix, ordre ligne par ligne) ; numero est une permutation de 0..M-1. */
typedef enum {
    NUMEROTATION_LIGNE = 0,      /* numero[i] = i */
    NUMEROTATION_RCMK,           /* Cuthill-McKee inverse */
    NUMEROTATION_DISSECTION,     /* dissection emboitee geometrique (grille) */
    NUMEROTATION_METIS           /* METIS_NodeND (compiler avec METIS=1) */
} MethodeNumerotation;

/* Nom lisible d'une methode (pour les affichages). */
const char *nomNumerotation(MethodeNumerotation methode);

/* Lit "ligne", "rcmk", "nd" ou "metis". Renvoie 0 en cas de succes, -1 sinon. */
int lireNumerotation(const char *texte, MethodeNumerotation *methode);

/* Calcule numero[] pour la grille de nxi x nyi noeuds interieurs.
   Renvoie 0 en cas de succes, -1 en cas d'echec (allocation, METIS absent...). */
int calculerNumerotation(MethodeNumerotation methode, int nxi, int nyi, int numero[]);

/* RCMK sur le graphe d'une matrice symetrique (triangle inferieur, CSC).
   numero[i] = nouvelle position de l'inconnue i. Renvoie 0 ou -1. */
int numerotationRCMK(const MatriceCSC *A, int numero[]);

#endif