#include <stdlib.h>
#include "cholesky.h"

static void libererListes(int **liste, size_t *taille, int M)
{
    if (liste != NULL)
        for (int j = 0; j < M; j++) free(liste[j]);
    free(liste);
    free(taille);
}

/*
 * A : triangle inferieur en CSC, lignes triees dans chaque colonne.
 * liste[j] = voisins de j d'indice > j (voisins de A + remplissage), tries, sans doublon.
 * En cas d'echec d'allocation : L->n = 0, L->nnz = 0, tableaux NULL.
 */
void choleskySymbolique(const MatriceCSC *A, MatriceCSC *L)
{
    L->n = 0;
    L->nnz = 0;
    L->valeurs = NULL;
    L->lignes = NULL;
    L->debutColonne = NULL;

    const int M = A->n;
    if (M <= 0) return;

    int **liste = calloc((size_t)M, sizeof(int *));
    size_t *taille = calloc((size_t)M, sizeof(size_t));
    if (liste == NULL || taille == NULL) {
        libererListes(liste, taille, M);
        return;
    }

    /* 1. Graphe initial : on ne garde que les voisins d'indice superieur. */
    for (int j = 0; j < M; j++) {
        size_t deb = A->debutColonne[j], fin = A->debutColonne[j + 1];
        size_t s = 0;
        for (size_t p = deb; p < fin; p++)
            if (A->lignes[p] > j) s++;
        if (s > 0) {
            liste[j] = malloc(s * sizeof(int));
            if (liste[j] == NULL) { libererListes(liste, taille, M); return; }
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
        if (s < 2) continue;

        int parent = liste[j][0];
        const int *S = liste[j] + 1;
        size_t reste = s - 1;
        size_t tp = taille[parent];

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

    /* 3. Passage au CSC : debutColonne[j+1] = debutColonne[j] + 1 + |S_j|. */
    size_t *deb = malloc(((size_t)M + 1) * sizeof(size_t));
    if (deb == NULL) { libererListes(liste, taille, M); return; }
    deb[0] = 0;
    for (int j = 0; j < M; j++)
        deb[j + 1] = deb[j] + 1 + taille[j];
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
        lignes[k++] = j;                         /* diagonale en premier */
        for (size_t t = 0; t < taille[j]; t++)
            lignes[k++] = liste[j][t];
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
 *
 * Colonne j :
 *   1. on depose la colonne j de A dans un vecteur de travail dense w (taille M) ;
 *   2. pour chaque colonne k < j qui contient une ligne j (l_jk != 0), on retire
 *      l_jk * (colonne k a partir de la ligne j) de w ;
 *   3. l_jj = sqrt(w[j]), puis l_ij = w[i] / l_jj pour les lignes i de la colonne j de L ;
 *   4. on remet w a zero sur la structure de la colonne j.
 *
 * Pour trouver les colonnes k de l'etape 2 sans parcourir les lignes de L (stockee par colonnes),
 * on utilise des listes chainees : tete[j] = premiere colonne k dont la prochaine ligne
 * a traiter est j, suivant[k] = colonne suivante dans la meme liste, pos[k] = position
 * courante dans la colonne k. Chaque coefficient de L est ainsi visite un nombre de fois
 * proportionnel a la taille de sa colonne, sans jamais construire L^T ni de matrice dense.
 *
 * Cela suppose que la structure de L est celle de l'analyse symbolique : pour une colonne k,
 * toutes les lignes i > j de la colonne k sont aussi dans la colonne j (propriete du remplissage).
 */
int choleskyNumerique(const MatriceCSC *A, MatriceCSC *L)
{
    if (A == NULL || L == NULL) return -1;
    const int M = A->n;
    if (M <= 0 || L->n != M || L->valeurs == NULL ||
        L->lignes == NULL || L->debutColonne == NULL) return -1;

    double *w      = calloc((size_t)M, sizeof(double));
    int    *tete   = malloc((size_t)M * sizeof(int));
    int    *suivant = malloc((size_t)M * sizeof(int));
    size_t *pos    = malloc((size_t)M * sizeof(size_t));
    if (w == NULL || tete == NULL || suivant == NULL || pos == NULL) {
        free(w); free(tete); free(suivant); free(pos);
        return -1;
    }
    for (int i = 0; i < M; i++) tete[i] = -1;

    int code = 0;
    for (int j = 0; j < M; j++) {
        const size_t deb = L->debutColonne[j];
        const size_t fin = L->debutColonne[j + 1];

        /* 1. Colonne j de A (triangle inferieur : lignes >= j). */
        for (size_t p = A->debutColonne[j]; p < A->debutColonne[j + 1]; p++)
            w[A->lignes[p]] = A->valeurs[p];

        /* 2. Contributions des colonnes k < j qui touchent la ligne j. */
        int k = tete[j];
        while (k != -1) {
            const int suiv = suivant[k];
            const size_t p = pos[k];                 /* ligne de la position p = j */
            const size_t finK = L->debutColonne[k + 1];
            const double ljk = L->valeurs[p];

            for (size_t q = p; q < finK; q++)
                w[L->lignes[q]] -= L->valeurs[q] * ljk;

            pos[k] = p + 1;
            if (p + 1 < finK) {                      /* k attend maintenant sa ligne suivante */
                const int r = L->lignes[p + 1];
                suivant[k] = tete[r];
                tete[r] = k;
            }
            k = suiv;
        }

        /* 3. Pivot puis colonne j de L. */
        const double d = w[j];
        if (!(d > 0.0)) {                            /* vrai aussi pour NaN */
            code = j + 1;
            break;
        }
        const double ljj = sqrt(d);
        L->valeurs[deb] = ljj;
        w[j] = 0.0;
        for (size_t p = deb + 1; p < fin; p++) {
            const int i = L->lignes[p];
            L->valeurs[p] = w[i] / ljj;
            w[i] = 0.0;                              /* 4. remise a zero */
        }

        /* La colonne j devra servir pour sa premiere ligne sous la diagonale. */
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

