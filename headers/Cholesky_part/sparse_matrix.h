#include <stddef.h>
typedef struct {
int n; /* Matrice carree de taille n. */
size_t nnz; /* Nombre de positions stockees. */
double *valeurs; /* nnz coefficients. */
int *lignes; /* nnz indices de ligne. */
size_t *debutColonne; /* n+1 positions dans les tableaux. */
} MatriceCSC;

/*Allouer une matrice sous la forme creuse*/
MatriceCSC *allocate_matrix(int n, size_t nnz);
/*Libérer la matrice creuse*/
void free_matrix(MatriceCSC *mat);

