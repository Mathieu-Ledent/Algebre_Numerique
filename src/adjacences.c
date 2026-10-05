#include "adjacences.h"
#include <stdlib.h>
#include <stdio.h>

#include <flann/flann.h>


void calculerAdjacences(const double points[][2], Adjacences *adj){
    //adj->N = nombre total de points
    int N = adj->N;


    //adj->k = nombre de voisins que l'on veut
    int k = adj->k;

    int nbResultats = k +1;

    struct FLANNParameters params = DEFAULT_FLANN_PARAMETERS;  //les réglages pour FLANN

    //On utilise le kd-tree :
    params.algorithm = FLANN_INDEX_KDTREE;   
    params.trees = 4;
    params.checks = FLANN_CHECKS_UNLIMITED;   //Correspond à une recherche lente mais très précise

    double *dataset = (double *)points;  //La raison de la conversion est car FLANN demande un double et on lui donne un const double

    float speedup;

    flann_index_t index = flann_build_index_double(dataset, N, 2, &speedup, &params);
    
    int *indicesResultats = malloc(nbResultats*sizeof(int));
    double *distances = malloc(nbResultats*sizeof(double));

    if(indicesResultats == NULL || distances == NULL){
        free(indicesResultats);
        free(distances);

        flann_free_index_double(index, &params);

        return;
    }

    for(int i=0; i<N; i++){
        double *query = dataset + 2 * i;  // -> c'est le point i auquel on va chercher ses voisins

        int status = flann_find_nearest_neighbors_index_double(index, query, 1, indicesResultats, distances, nbResultats, &params);

        if(status != 0){
            fprintf(stderr, "Erreur FLANN pour le point %d.\n", i);
            continue;
        }

        int compteur = 0;       // On s'assure que le point dont l'index est i n'est pas pris en compte

        for(int j =0; j<nbResultats && compteur<k; j++){
            int voisin = indicesResultats[j];

            if(voisin != i){
                adj->voisins[i*k+compteur] = voisin;

                compteur++;
            }
        }
    }

    flann_free_index_double(index, &params);

    free(indicesResultats);
    free(distances);


}