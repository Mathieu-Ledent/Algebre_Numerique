# LINMA1170 - Devoir 1 : Factorisation LU et Cholesky

**Année académique :** 2025-2026  
**Équipe enseignante :** V. Degrooff, C. Heneffe, N. Tihon, J.-F. Remacle

Ce dépôt contient les implémentations et les analyses numériques requises pour le premier devoir du cours LINMA1170. Le projet se concentre sur l'implémentation bas niveau en C de deux décompositions fondamentales de l'algèbre linéaire : la factorisation LU avec pivotage complet et la factorisation *fast* Cholesky.

## 📁 Structure du projet

- `devoir_1.c` : Fichier source principal contenant les implémentations des algorithmes de factorisation en C.
- `analyses.py`  : Script générant les graphiques de complexité temporelle et l'analyse du *growth factor*.
- `presentation.pdf` : Slides de la présentation des résultats.

## 🚀 Implémentations (C)

Les matrices sont traitées sous forme de vecteurs unidimensionnels (*row-major*).

### 1. Factorisation LU avec pivotage complet
Algorithme de factorisation d'une matrice $A$ de taille $n \times n$. Contrairement au pivotage partiel, le pivot à l'étape $k$ est la plus grande entrée de la sous-matrice $A_{k:n,k:n}$, ce qui maximise la stabilité. 

Mathématiquement, on applique des permutations à gauche ($P$) et à droite ($Q$) :
$$PAQ = LU$$

**Signature :**
```c
int lu_factorisation(double *A, int *P, int *Q, int n);
```

### 2. Factorisation "Fast Cholesky"
Factorisation pour les matrices symétriques définies positives (SDP). Cette implémentation exploite l'algorithme de Strassen pour les produits matriciels afin d'accélérer les calculs, en ne travaillant idéalement que sur la partie triangulaire de la matrice.

**Signature :**
```c
int fast_Cholesky(double *A, int n);
```

## 📊 Analyses Numériques (Python/Script)

Le script de test associé réalise les expériences suivantes pour l'évaluation :

1. **Complexité temporelle** : Tracé d'un graphique en échelle *log-log* comparant les temps d'exécution de nos implémentations en fonction de la taille $n$, avec une comparaison par rapport aux librairies de référence (SciPy ou LAPACK).
2. **Growth Factor (Stabilité)** : Analyse du facteur de croissance $\rho$ sur des matrices aléatoires, comparant les stratégies sans pivotage, avec pivotage partiel, et avec pivotage complet :
   $$\rho = \frac{\max_{i,j}|U_{ij}|}{\max_{i,j}|A_{ij}|}$$
3. **(Bonus) Phénomène de Fill-in** : Analyse de l'augmentation des entrées non-nulles (sans pivotage) sur la matrice de l'opérateur Laplacien 2D discrétisé, en illustrant l'impact de la numérotation des nœuds.

## 🛠️ Compilation

Le code peut utiliser la librairie externe BLAS (optionnel mais recommandé). Exemple de compilation type avec GCC :

```bash
gcc -Wall -O3 devoir_1.c -o devoir_1 -lblas -lm
```

## 📅 Échéances (Soumission via Gradescope)

- ⏳ **16 mars 2026 à 14h00 :** Soumission du code (`devoir_1.c`) et du script de reproduction des figures.
- ⏳ **23 mars 2026 à 12h00 :** Soumission des slides de présentation (`presentation.pdf`).

## 👥 Auteurs

- **[Ledent Mathieu]**
