#include <stdlib.h>
#include <string.h>
#include "numerotation.h"
#include "poisson.h"

#ifdef AVEC_METIS
#include <metis.h>
#endif

const char *nomNumerotation(MethodeNumerotation methode)
{
    switch (methode) {
    case NUMEROTATION_LIGNE:      return "ligne par ligne";
    case NUMEROTATION_RCMK:       return "RCMK";
    case NUMEROTATION_DISSECTION: return "dissection emboitee";
    case NUMEROTATION_METIS:      return "METIS";
    }
    return "?";
}

int lireNumerotation(const char *texte, MethodeNumerotation *methode)
{
    if (strcmp(texte, "ligne") == 0)      *methode = NUMEROTATION_LIGNE;
    else if (strcmp(texte, "rcmk") == 0)  *methode = NUMEROTATION_RCMK;
    else if (strcmp(texte, "nd") == 0)    *methode = NUMEROTATION_DISSECTION;
    else if (strcmp(texte, "metis") == 0) *methode = NUMEROTATION_METIS;
    else return -1;
    return 0;
}

/* Graphe symetrique (CSR, sans la diagonale) a partir du triangle inferieur de A. */
static int grapheSymetrique(const MatriceCSC *A, int **debutOut, int **voisinsOut)
{
    const int M = A->n;
    int *deb = calloc((size_t)M + 1, sizeof(int));
    if (deb == NULL) return -1;
    for (int j = 0; j < M; j++)
        for (size_t p = A->debutColonne[j]; p < A->debutColonne[j + 1]; p++) {
            int i = A->lignes[p];
            if (i != j) { deb[i + 1]++; deb[j + 1]++; }
        }
    for (int i = 0; i < M; i++) deb[i + 1] += deb[i];

    int *vois = malloc(((size_t)deb[M] > 0 ? (size_t)deb[M] : 1) * sizeof(int));
    int *rempli = malloc((size_t)M * sizeof(int));
    if (vois == NULL || rempli == NULL) { free(deb); free(vois); free(rempli); return -1; }
    memcpy(rempli, deb, (size_t)M * sizeof(int));
    for (int j = 0; j < M; j++)
        for (size_t p = A->debutColonne[j]; p < A->debutColonne[j + 1]; p++) {
            int i = A->lignes[p];
            if (i != j) { vois[rempli[i]++] = j; vois[rempli[j]++] = i; }
        }
    free(rempli);
    *debutOut = deb;
    *voisinsOut = vois;
    return 0;
}

/* Parcours en largeur depuis 'depart' (sommets non marques uniquement), voisins par degre
   croissant. Ecrit l'ordre de visite dans ordre[] a partir de n0 et renvoie le nouveau total.
   niveau[] recoit la distance au depart (utile pour le noeud pseudo-peripherique). */
static int parcoursCM(const int *deb, const int *vois, int depart, int *marque, int marqueur,
                      int *ordre, int n0, int *niveau)
{
    int tete = n0, queue = n0;
    ordre[queue++] = depart;
    marque[depart] = marqueur;
    if (niveau) niveau[depart] = 0;
    while (tete < queue) {
        int i = ordre[tete++];
        int debutAjout = queue;
        for (int p = deb[i]; p < deb[i + 1]; p++) {
            int j = vois[p];
            if (marque[j] != marqueur) {
                marque[j] = marqueur;
                if (niveau) niveau[j] = niveau[i] + 1;
                ordre[queue++] = j;
            }
        }
        /* Tri par insertion des voisins ajoutes, par degre croissant (au plus quelques uns). */
        for (int a = debutAjout + 1; a < queue; a++) {
            int x = ordre[a], dx = deb[x + 1] - deb[x], b = a;
            while (b > debutAjout && deb[ordre[b - 1] + 1] - deb[ordre[b - 1]] > dx) {
                ordre[b] = ordre[b - 1];
                b--;
            }
            ordre[b] = x;
        }
    }
    return queue;
}

int numerotationRCMK(const MatriceCSC *A, int numero[])
{
    const int M = A->n;
    if (M <= 0) return -1;
    int *deb, *vois;
    if (grapheSymetrique(A, &deb, &vois) != 0) return -1;

    int *ordre  = malloc((size_t)M * sizeof(int));
    int *marque = calloc((size_t)M, sizeof(int));   /* 0 : libre, 1 : numerote, >= 2 : temporaire */
    int *niveau = malloc((size_t)M * sizeof(int));
    int *tmp    = malloc((size_t)M * sizeof(int));
    if (!ordre || !marque || !niveau || !tmp) {
        free(deb); free(vois); free(ordre); free(marque); free(niveau); free(tmp);
        return -1;
    }

    int total = 0, marqueur = 2;
    for (int s = 0; s < M; s++) {
        if (marque[s] == 1) continue;

        /* Noeud de depart de degre minimal dans la composante de s... */
        int nComp = parcoursCM(deb, vois, s, marque, marqueur++, tmp, 0, NULL);
        int depart = s;
        for (int a = 0; a < nComp; a++) {
            int x = tmp[a];
            if (deb[x + 1] - deb[x] < deb[depart + 1] - deb[depart]) depart = x;
        }
        /* ...puis recherche d'un noeud pseudo-peripherique (George-Liu) : on repart du noeud
           de degre minimal du dernier niveau tant que l'excentricite augmente. */
        int excentricite = -1;
        for (int essai = 0; essai < 10; essai++) {
            parcoursCM(deb, vois, depart, marque, marqueur++, tmp, 0, niveau);
            int dernier = tmp[nComp - 1], e = niveau[dernier];
            if (e <= excentricite) break;
            excentricite = e;
            int meilleur = dernier;
            for (int a = nComp - 1; a >= 0 && niveau[tmp[a]] == e; a--)
                if (deb[tmp[a] + 1] - deb[tmp[a]] < deb[meilleur + 1] - deb[meilleur]) meilleur = tmp[a];
            depart = meilleur;
        }

        /* Cuthill-McKee definitif depuis le noeud retenu. */
        total = parcoursCM(deb, vois, depart, marque, 1, ordre, total, NULL);
    }

    /* Inversion de l'ordre : ordre[t] est la t-ieme inconnue de CM, elle recoit M-1-t. */
    for (int t = 0; t < M; t++) numero[ordre[t]] = M - 1 - t;

    free(deb); free(vois); free(ordre); free(marque); free(niveau); free(tmp);
    return 0;
}

/* Dissection emboitee geometrique sur le rectangle [x0,x1) x [y0,y1) de noeuds interieurs :
   on numerote recursivement les deux moities, puis la ligne separatrice en dernier. */
static void dissection(int nxi, int x0, int x1, int y0, int y1, int numero[], int *compteur)
{
    const int lx = x1 - x0, ly = y1 - y0;
    if (lx <= 0 || ly <= 0) return;
    if (lx <= 2 || ly <= 2) {
        for (int iy = y0; iy < y1; iy++)
            for (int ix = x0; ix < x1; ix++)
                numero[iy * nxi + ix] = (*compteur)++;
        return;
    }
    if (lx >= ly) {
        int m = x0 + lx / 2;
        dissection(nxi, x0, m, y0, y1, numero, compteur);
        dissection(nxi, m + 1, x1, y0, y1, numero, compteur);
        for (int iy = y0; iy < y1; iy++) numero[iy * nxi + m] = (*compteur)++;
    } else {
        int m = y0 + ly / 2;
        dissection(nxi, x0, x1, y0, m, numero, compteur);
        dissection(nxi, x0, x1, m + 1, y1, numero, compteur);
        for (int ix = x0; ix < x1; ix++) numero[m * nxi + ix] = (*compteur)++;
    }
}

#ifdef AVEC_METIS
static int numerotationMETIS(const MatriceCSC *A, int numero[])
{
    const int M = A->n;
    int *deb, *vois;
    if (grapheSymetrique(A, &deb, &vois) != 0) return -1;
    idx_t n = M;
    idx_t *xadj = malloc(((size_t)M + 1) * sizeof(idx_t));
    idx_t *adjncy = malloc(((size_t)deb[M] > 0 ? (size_t)deb[M] : 1) * sizeof(idx_t));
    idx_t *perm = malloc((size_t)M * sizeof(idx_t));
    idx_t *iperm = malloc((size_t)M * sizeof(idx_t));
    int code = -1;
    if (xadj && adjncy && perm && iperm) {
        for (int i = 0; i <= M; i++) xadj[i] = deb[i];
        for (int p = 0; p < deb[M]; p++) adjncy[p] = vois[p];
        idx_t options[METIS_NOPTIONS];
        METIS_SetDefaultOptions(options);
        options[METIS_OPTION_NUMBERING] = 0;
        if (METIS_NodeND(&n, xadj, adjncy, NULL, options, perm, iperm) == METIS_OK) {
            /* iperm[i] = position de l'inconnue i dans la matrice permutee. */
            for (int i = 0; i < M; i++) numero[i] = (int)iperm[i];
            code = 0;
        }
    }
    free(deb); free(vois); free(xadj); free(adjncy); free(perm); free(iperm);
    return code;
}
#endif

int calculerNumerotation(MethodeNumerotation methode, int nxi, int nyi, int numero[])
{
    if (nxi <= 0 || nyi <= 0 || numero == NULL) return -1;
    const int M = nxi * nyi;

    if (methode == NUMEROTATION_LIGNE) {
        for (int i = 0; i < M; i++) numero[i] = i;
        return 0;
    }
    if (methode == NUMEROTATION_DISSECTION) {
        int compteur = 0;
        dissection(nxi, 0, nxi, 0, nyi, numero, &compteur);
        return 0;
    }

    /* RCMK et METIS travaillent sur le graphe de la matrice assemblee ligne par ligne. */
    for (int i = 0; i < M; i++) numero[i] = i;
    Grille g = { nxi + 2, nyi + 2, 0.0, 0.0, 1.0, NULL, NULL };
    MatriceCSC A = {0};
    construireMatricePoisson(&g, numero, &A);
    if (A.n != M) return -1;

    int code = -1;
    if (methode == NUMEROTATION_RCMK)
        code = numerotationRCMK(&A, numero);
#ifdef AVEC_METIS
    else if (methode == NUMEROTATION_METIS)
        code = numerotationMETIS(&A, numero);
#endif
    libererTableauxCSC(&A);
    return code;
}
