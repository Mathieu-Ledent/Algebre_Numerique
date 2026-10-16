#include <stdio.h>
#include <stdlib.h>
#include "poisson.h"
#include "cholesky.h"

/* Compare la structure de L avec un remplissage de reference calcule par
   elimination booleenne dense (independante de choleskySymbolique). */
static int verifierStructure(const MatriceCSC *A, const MatriceCSC *L)
{
    const int M = A->n;
    int ok = 1;
    if (L->n != M) { printf("  ECHEC : L->n incorrect\n"); return 0; }
    if (L->debutColonne[0] != 0 || L->debutColonne[M] != L->nnz) {
        printf("  ECHEC : debutColonne incorrect\n"); return 0;
    }

    char *fill = calloc((size_t)M * M, 1);
    for (int j = 0; j < M; j++) fill[(size_t)j * M + j] = 1;
    for (int c = 0; c < M; c++)
        for (size_t p = A->debutColonne[c]; p < A->debutColonne[c + 1]; p++)
            fill[(size_t)A->lignes[p] * M + c] = 1;

    for (int j = 0; j < M; j++)
        for (int u = j + 1; u < M; u++)
            if (fill[(size_t)u * M + j])
                for (int v = u + 1; v < M; v++)
                    if (fill[(size_t)v * M + j]) fill[(size_t)v * M + u] = 1;

    size_t attendu = 0;
    for (int r = 0; r < M; r++)
        for (int c = 0; c <= r; c++) attendu += fill[(size_t)r * M + c];

    if (attendu != L->nnz) {
        printf("  ECHEC : nnz(L) = %zu, attendu %zu\n", L->nnz, attendu);
        ok = 0;
    }
    for (int c = 0; c < M && ok; c++) {
        size_t deb = L->debutColonne[c], fin = L->debutColonne[c + 1];
        if (fin <= deb || L->lignes[deb] != c) {
            printf("  ECHEC : diagonale absente en tete de la colonne %d\n", c); ok = 0; break;
        }
        for (size_t p = deb; p < fin; p++) {
            int r = L->lignes[p];
            if (r < c || r >= M || !fill[(size_t)r * M + c]) {
                printf("  ECHEC : position (%d,%d) inattendue\n", r, c); ok = 0; break;
            }
            if (p > deb && L->lignes[p] <= L->lignes[p - 1]) {
                printf("  ECHEC : lignes non croissantes (colonne %d)\n", c); ok = 0; break;
            }
        }
    }
    free(fill);
    return ok;
}

/* mode 0 : identite ; 1 : inversee ; 2 : (7 i) % M (permutation si pgcd(7, M) = 1).
   verifierFormule : pour la numerotation ligne par ligne, la bande est entierement remplie. */
static int tester(const char *nom, int nx, int ny, int mode, int verifierFormule)
{
    Grille g = { nx, ny, 0.0, 0.0, 1.0, NULL, NULL };
    const int nxi = nx - 2, M = nxi * (ny - 2);
    int *numero = malloc((size_t)M * sizeof(int));
    for (int i = 0; i < M; i++)
        numero[i] = (mode == 0) ? i : (mode == 1) ? M - 1 - i : (7 * i) % M;

    MatriceCSC A = {0}, L = {0};
    construireMatricePoisson(&g, numero, &A);
    choleskySymbolique(&A, &L);

    int ok = verifierStructure(&A, &L);

        if (ok && verifierFormule && (ny - 2) >= 2) {   /* formule valable pour au moins 2 rangees */
        size_t formule = 0;
        for (int j = 0; j < M; j++) {
            int sous = (M - 1 - j < nxi) ? M - 1 - j : nxi;
            formule += 1 + (size_t)sous;
        }
        if (nxi >= 2)
            formule -= (size_t)(nxi - 1) * (size_t)(nxi - 2) / 2;   /* premiere rangee : pas de bande pleine */
        if (formule != L.nnz) {
            printf("  ECHEC : formule de bande = %zu, obtenu %zu\n", formule, L.nnz);
            ok = 0;
        }
    }
    printf("%s : M = %d, nnz(A) = %zu, nnz(L) = %zu -> %s\n",
           nom, M, A.nnz, L.nnz, ok ? "OK" : "ECHEC");

    libererTableauxCSC(&A);
    libererTableauxCSC(&L);
    free(numero);
    return ok;
}

int main(void)
{
    /* Exemple de l'enonce : aretes (0,1), (0,2), (1,3), (2,3), resultat attendu : 9 positions. */
    int tout = 1;
    {
        double v[] = {4, -1, -1, 4, -1, 4, -1, 4};
        int    l[] = {0, 1, 2, 1, 3, 2, 3, 3};
        size_t d[] = {0, 3, 5, 7, 8};
        MatriceCSC A = { 4, 8, v, l, d };
        MatriceCSC L = {0};
        choleskySymbolique(&A, &L);

        size_t dAtt[] = {0, 3, 6, 8, 9};
        int    lAtt[] = {0, 1, 2, 1, 2, 3, 2, 3, 3};
        int ok = (L.nnz == 9);
        for (int j = 0; j <= 4 && ok; j++) ok = (L.debutColonne[j] == dAtt[j]);
        for (size_t k = 0; k < 9 && ok; k++) ok = (L.lignes[k] == lAtt[k]);
        printf("Exemple de l'enonce (4 inconnues) : %s\n", ok ? "OK" : "ECHEC");
        tout &= ok;
        libererTableauxCSC(&L);   /* A pointe vers des tableaux locaux : on ne la libere pas */
    }

    tout &= tester("4x4   ligne par ligne ", 4, 4, 0, 1);
    tout &= tester("8x6   ligne par ligne ", 8, 6, 0, 1);
    tout &= tester("22x22 ligne par ligne ", 22, 22, 0, 1);
    tout &= tester("8x6   inversee        ", 8, 6, 1, 0);
    tout &= tester("12x12 permutation 7i  ", 12, 12, 2, 0);   /* M = 100, pgcd(7, 100) = 1 */
    tout &= tester("3x3   (1 inconnue)    ", 3, 3, 0, 1);

    printf("\n%s\n", tout ? "Tous les tests passent." : "Au moins un test echoue.");
    return tout ? 0 : 1;
}