#ifndef ADJACENCES_H
#define ADJACENCES_H

typedef struct{
    int N;
    int k;
    int *voisins;
} Adjacences;

void calculerAdjacences(const double points[][2], Adjacences *adj);

#endif