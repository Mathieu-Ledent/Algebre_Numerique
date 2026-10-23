#include <math.h>
#include <stdlib.h>
#include "cholesky.h"

/* Optimisations GCC des noyaux de calcul (ignorees par les autres compilateurs).
   Pas de -ffast-math : on garde l'arithmetique IEEE stricte pour Cholesky. */
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize("O3", "unroll-loops")
#endif

static void libererListes(int **liste, size_t *taille, int M)
{
    if (liste != NULL)
        for (int j = 0; j < M; j++) free(liste[j]);
    free(liste);
    free(taille);
}

/*
 * A : matrice stockée en CSC, L contient les positions à remplir dans le triangle inférieur 
 */
void choleskySymbolique(const MatriceCSC *A, MatriceCSC *L)
{
    /* Initialisation de la matrice de sortie L*/
    L->n = 0;
    L->nnz = 0;
    L->valeurs = NULL;
    L->lignes = NULL;
    L->debutColonne = NULL;

    const int M = A->n;
    if (M <= 0) return;

    /* Chaque noeud possède une liste de voisins non encore*/
    int **liste = calloc((size_t)M, sizeof(int *));
    size_t *taille = calloc((size_t)M, sizeof(size_t));
    if (liste == NULL || taille == NULL) {
        libererListes(liste, taille, M);
        return;
    }

    /* On construit le graphe initial à partir de A */
    for (int j = 0; j < M; j++) {
        size_t deb = A->debutColonne[j], fin = A->debutColonne[j + 1];
        size_t s = 0; /* Nombre de voisins*/
        for (size_t p = deb; p < fin; p++)
            if (A->lignes[p] > j) s++; /* Eléments hors diagonale ont un indice supérieur à j*/
        if (s > 0) {
            liste[j] = malloc(s * sizeof(int));
            if (liste[j] == NULL) { libererListes(liste, taille, M); return; }
            /*On remplit la liste avec les voisins de j*/
            size_t k = 0;
            for (size_t p = deb; p < fin; p++)
                if (A->lignes[p] > j) liste[j][k++] = A->lignes[p];
            taille[j] = s;
        }
    }

    /* 2. Elimination symbolique.
       Quand j est elimine, ses voisins restants S_j doivent tous etre relies entre eux.
       Il suffit de fusionner S_j \ {min} dans la liste de min (le "parent" de j) :
       les autres liaisons seront propagees quand le parent sera elimine a son tour. */
    for (int j = 0; j < M; j++) {
        size_t s = taille[j];
        if (s < 2) continue; /*impossible de créer une nouvelle connexion avec moins de 2 voisins */

        int parent = liste[j][0]; /* premier voisin de la liste*/
        const int *S = liste[j] + 1;
        size_t reste = s - 1; /* nombre de voisins restants*/
        size_t tp = taille[parent]; /* taille de la liste actuelle du parent*/

        int *fusion = malloc((tp + reste) * sizeof(int));
        if (fusion == NULL) { libererListes(liste, taille, M); return; }

        size_t a = 0, b = 0, k = 0;
        while (a < tp && b < reste) {
            int x = liste[parent][a], y = S[b];
            if (x < y)      { fusion[k++] = x; a++; }
            else if (x > y) { fusion[k++] = y; b++; }
            else            { fusion[k++] = x; a++; b++; }
        }
        while (a < tp)    fusion[k++] = liste[parent][a++];
        while (b < reste) fusion[k++] = S[b++];

        free(liste[parent]);
        liste[parent] = fusion;
        taille[parent] = k;
    }

    /* On stocke ces positions en format CSC : debutColonne[j+1] = debutColonne[j] + 1 + |S_j|. */
    size_t *deb = malloc(((size_t)M + 1) * sizeof(size_t));
    if (deb == NULL) { libererListes(liste, taille, M); return; }
    deb[0] = 0;
    for (int j = 0; j < M; j++)
        deb[j + 1] = deb[j] + 1 + taille[j]; /*Chaque colonne de L possède un élément diagonal et taille[j] éléments hors diagonale*/
    const size_t nnz = deb[M];

    int *lignes = malloc(nnz * sizeof(int));
    double *valeurs = calloc(nnz, sizeof(double));
    if (lignes == NULL || valeurs == NULL) {
        free(deb); free(lignes); free(valeurs);
        libererListes(liste, taille, M);
        return;
    }

    for (int j = 0; j < M; j++) {
        size_t k = deb[j];
        lignes[k++] = j;                         /* on remplit la diagonale en premier */
        for (size_t t = 0; t < taille[j]; t++)
            lignes[k++] = liste[j][t]; /* puis les éléments hors diagonale*/
    }

    libererListes(liste, taille, M);

    L->n = M;
    L->nnz = nnz;
    L->debutColonne = deb;
    L->lignes = lignes;
    L->valeurs = valeurs;
}


/*
 * Cholesky "left-looking" creux, colonne par colonne.
 * Depuis les positions trouvées par la version symbolique on va calculer les coefficients.
 */
int choleskyNumerique(const MatriceCSC *A, MatriceCSC *L)
{
    if (A == NULL || L == NULL) return -1;
    const int M = A->n;
    if (M <= 0 || L->n != M || L->valeurs == NULL ||
        L->lignes == NULL || L->debutColonne == NULL) return -1;

    double *w      = calloc((size_t)M, sizeof(double)); /*contient les coeff mis à jour pour la colonne en cours*/
    int    *tete   = malloc((size_t)M * sizeof(int)); /*1ère colonne de L à traiter pour une ligne donnée*/
    int    *suivant = malloc((size_t)M * sizeof(int)); 
    size_t *pos    = malloc((size_t)M * sizeof(size_t));/*Pour la colonne active, donne la position du prochain coeff à traiter*/
    if (w == NULL || tete == NULL || suivant == NULL || pos == NULL) {
        free(w); free(tete); free(suivant); free(pos);
        return -1;
    }
    /*On travaille avec des listes chaînées*/
    for (int i = 0; i < M; i++) tete[i] = -1;

    int code = 0;
    /*On parcourt la colonne j*/
    for (int j = 0; j < M; j++) {
        const size_t deb = L->debutColonne[j];
        const size_t fin = L->debutColonne[j + 1];

        /* On copie la colonne j de A dans w */
        for (size_t p = A->debutColonne[j]; p < A->debutColonne[j + 1]; p++)
            w[A->lignes[p]] = A->valeurs[p];

        /* Contributions des colonnes précédentes */
        int k = tete[j];
        while (k != -1) {
            const int suiv = suivant[k];
            const size_t p = pos[k];                 /* ligne de la position p = j */
            const size_t finK = L->debutColonne[k + 1];
            const double ljk = L->valeurs[p];

            for (size_t q = p; q < finK; q++)
            /*Pour chaque ligne i à partir de j, on effectue la mise à jour*/
                w[L->lignes[q]] -= L->valeurs[q] * ljk;

            pos[k] = p + 1;
            if (p + 1 < finK) {                      /* k attend maintenant sa ligne suivante */
                const int r = L->lignes[p + 1];
                suivant[k] = tete[r];
                tete[r] = k;
            }
            k = suiv;
        }

        /* Calculs le pivot diagonal. */
        const double d = w[j];
        if (!(d > 0.0)) {                          
            code = j + 1;
            break;
        }
        const double ljj = sqrt(d);
        L->valeurs[deb] = ljj;
        w[j] = 0.0;/* Remise à zéro pour les calculs suivants*/
        for (size_t p = deb + 1; p < fin; p++) { /*Calcul des éléments sous la diagonale*/
            const int i = L->lignes[p]; 
            L->valeurs[p] = w[i] / ljj;
            w[i] = 0.0;                              /* Remise à zéro pour les calculs suivants*/
        }

        /* La colonne j devra contribuer au calcul d'une ligne future (si elle possède >1 élément sous la diagonale) */
        if (fin - deb > 1) {
            pos[j] = deb + 1;
            const int r = L->lignes[deb + 1];
            suivant[j] = tete[r];
            tete[r] = j;
        }
    }

    free(w); free(tete); free(suivant); free(pos);
    return code;
}

