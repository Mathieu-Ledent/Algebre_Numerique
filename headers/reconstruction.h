#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H

#include <stddef.h>
#include "grille.h"
#include "numerotation.h"

/* Options facultatives, a regler avant l'appel (valeurs par defaut : ligne par ligne, a_i = 1). */
typedef struct {
    MethodeNumerotation numerotation;
    int poidsAdaptes;          /* 0 : a_i = 1 ; 1 : poids selon l'espacement local */
    int inverserNormales;      /* 0 : rien ; 1 : toutes ; 2 : seulement la composante x (etude 6.2) */
} OptionsReconstruction;

/* Mesures de la derniere reconstruction (temps en secondes). */
typedef struct {
    double tAdjacences, tNormales, tOrientation, tDivergence;
    double tNumerotation, tAssemblage, tSymbolique, tNumerique, tSubstitutions, tNiveau, tTotal;
    int M;                     /* nombre d'inconnues */
    size_t nnzA, nnzL;         /* positions stockees (triangle inferieur) */
    int factoriseReutilise;    /* 1 si numerotation + facteur viennent du cache */
    double residuRelatif;      /* ||AF - b|| / ||b|| (ou residu absolu si b = 0) */
    double ecartTypeNiveau;    /* dispersion des f_h(p_i) autour de tau */
} StatistiquesReconstruction;

void reglerOptionsReconstruction(const OptionsReconstruction *options);
const StatistiquesReconstruction *statistiquesReconstruction(void);

/* Libere la numerotation et le facteur gardes en memoire entre deux appels. */
void libererCacheReconstruction(void);

/* Cahier des charges : ecrit f (nx*ny valeurs, f[r + nx*s], bord nul) et *niveau = tau.
   Renvoie 0 en cas de succes, une valeur non nulle sinon. */
int reconstructionPoisson(int N, const double points[][2],
                          int k, double sigma,
                          double xmin, double ymin, double h,
                          int nx, int ny,
                          double *f, double *niveau);

/* Exporte le nuage, la grille, le niveau et f dans un fichier texte :
     # nuage            puis N, puis N lignes "x y"
     # grille           puis "nx ny xmin ymin h"
     # niveau           puis tau
     # f                puis ny lignes de nx valeurs (ligne s : f[r + nx*s], r = 0..nx-1)
   Renvoie 0 en cas de succes. */
int exporterResultats(const char *fichier, int N,
                      const double points[][2],
                      const Grille *grille, double niveau);

#endif
