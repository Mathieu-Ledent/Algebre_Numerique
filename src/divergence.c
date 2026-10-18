#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct {
    int nx, ny; /* Dimensions, bord compris. */
    double xmin, ymin; /* Origine de la grille. */
    double h; /* Pas identique en x et y. */
    double *g; /* Divergence : nx*ny valeurs. */
    double *f; /* Solution : nx*ny valeurs. */
} Grille;

double guassianSplitting(double dist2, double sigma) {
    double coef = 1/(2*M_PI*sigma*sigma);
    double expo = exp((-1) * (dist2/(2*sigma*sigma)));
    return coef*expo;
}

void calculerDivergence(int N, const double points[][2],
    const double normales[][2],
    const double poids[], double sigma,
    Grille *grille) 
    {

    for (int j = 0; j < grille->ny; j++) {
        for (int i = 0; i < grille->nx; i++) {

            double div = 0;

            double qx = grille->xmin + i*grille->h;
            double qy = grille->ymin + j*grille->h;
        
            for (int k = 0; k < N; k++) {
                
                double diff_x = qx - points[k][0];
                double diff_y = qy - points[k][1];
                double dist2 = (diff_x*diff_x) + (diff_y*diff_y);

                double g_value = guassianSplitting(dist2, sigma);
                double sca_product = normales[k][0] * diff_x + normales[k][1] * diff_y;

                // Sécurité si le tableau poids est NULL (poids par défaut = 1)
                double p_k = (poids!=NULL) ? poids[k]:1;

                div += p_k * sca_product * g_value;
            }

            div *= (-1)/(sigma*sigma);
            int index = i + grille->nx*j;
            grille->g[index] = div;

        }
    }

    }