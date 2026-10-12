#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "adjacences.h"
#include "utils.h"

void CalculerDroitesMoyennes(const double points[][2], const Adjacences *adj,double barycentres[][2], double normales[][2]){
    int N = adj->N;
    int k = adj->k;

    for(int i=0; i<N; i++){
        double sommeX = points[i][0];
        double sommeY = points[i][1];

        for(int j=0; j<k; j++){
            int voisin = adj->voisins[i*k+j];

            sommeX += points[voisin][0];
            sommeY += points[voisin][1];
        }

        barycentres[i][0] = sommeX / (k+1);  //xi barre
        barycentres[i][1] = sommeY / (k+1);   //yi barre

        double Cxx = 0.0;
        double Cxy = 0.0;
        double Cyy = 0.0;

        double dx = points[i][0] - barycentres[i][0];
        double dy = points[i][1] - barycentres[i][1];

        Cxx += dx*dx;
        Cxy += dx*dy;
        Cyy += dy*dy;

        for(int j=0; j<k; j++){
            int voisin = adj->voisins[i*k +j];

            dx = points[voisin][0] - barycentres[voisin][0];
            dy = points[voisin][1] - barycentres[voisin][1];

            Cxx += dx*dx;
            Cxy += dx*dy;
            Cyy += dy*dy;
        }

        Cxx /= k+1;
        Cxy /= k+1;
        Cyy /= k+1;

        double delta = sqrt((Cxx-Cyy) * (Cxx-Cyy) + 4.0 * Cxy*Cxy);

        double lambdaMin = 0.5* (Cxx+Cyy-delta);

        double nx;
        double ny;

        if(fabs(Cxy) > 1e-12){
            nx = Cxy;
            ny = lambdaMin-Cxx;
        }else{
            if(Cxx <= Cyy){
                nx = 1.0;
                ny = 0.0;
            }else{
                nx = 0.0;
                ny = 1.0;
            }
        }

        double norme = sqrt(nx*nx + ny*ny);

        nx /= norme;
        ny /= norme;

        normales[i][0] = nx;
        normales[i][1] = ny;
    }
   
}



void orienterNormales(const Adjacences *adj, double normales[][2]) {
    int N = adj->N;
    int k = adj->k;
    
    // Tableau pour marquer les points "découverts" (0 = non, 1 = oui)
    int decouverts[N];
    memset(decouverts, 0, N * sizeof(int));
    
    // file (FIFO) pour stocker les indices des points à traiter
    int file[N]; 
    
    // On parcourt tous les points de 0 à N-1. 
    // Cela garantit que chaque composante connexe séparée sera traitée.
    for (int start_node = 0; start_node < N; start_node++) {
        
        // Si le point n'a pas encore été découvert, il devient un nouveau point de départ
        if (decouverts[start_node] == 0) {
            
            // Initialisation des pointeurs de la file d'attente
            int tete = 0;
            int queue = 0;
            
            // On ajoute le point de départ dans la file et on le marque
            file[queue++] = start_node;
            decouverts[start_node] = 1; 
            
            // Tant que la file n'est pas vide
            while (tete < queue) {
                // On "retire" le premier élément de la file pour le traiter
                int i = file[tete++];
                
                // On isole la normale du point courant i pour la comparaison
                double ref_x = normales[i][0];
                double ref_y = normales[i][1];
                
                // On examine tous ses k voisins
                for (int j = 0; j < k; j++) {
                    int indice_voisin = adj->voisins[i * k + j];
                    
                    // On n'agit que si ce voisin est "non encore découvert"
                    if (decouverts[indice_voisin] == 0) {
                        
                        // Calcul du produit scalaire 
                        double scal = (ref_x * normales[indice_voisin][0]) + 
                                      (ref_y * normales[indice_voisin][1]);
                        
                        // Si le produit scalaire est négatif, on retourne la flèche
                        if (scal < 0) {
                            normales[indice_voisin][0] = -normales[indice_voisin][0];
                            normales[indice_voisin][1] = -normales[indice_voisin][1];
                        }
                        
                        // On marque immédiatement ce voisin comme découvert
                        decouverts[indice_voisin] = 1;
                        
                        // Et on l'ajoute à l'arrière de la file pour qu'il propage son orientation plus tard
                        file[queue++] = indice_voisin;
                    }
                }
            }
        }
    }
}