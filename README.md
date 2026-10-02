# LINMA1170 - Devoir 1 : Reconstruction de surface (Poisson 2D)

**Année académique :** 2026-2027  
**Cours :** LINMA1170 - Analyse numérique  

Ce projet implémente une chaîne complète de reconstruction de surface (courbe 2D) à partir d'un nuage de points potentiellement bruité. La méthode repose sur la **reconstruction de Poisson**, combinée à un solveur de système linéaire utilisant une **factorisation de Cholesky creuse**.

## 🎯 Objectifs du projet

L'objectif est d'estimer une fonction implicite $f(x,y)$ sur une grille cartésienne régulière dont une ligne de niveau $\Gamma$ correspond à la courbe à reconstruire. Le projet est divisé en plusieurs modules qui s'enchaînent :

1. **Recherche spatiale (FLANN) :** Recherche des $k$ plus proches voisins pour chaque point du nuage.
2. **Estimation des normales :** Ajustement d'une droite moyenne par analyse de covariance spatiale locale pour déterminer la direction normale.
3. **Orientation cohérente (Graphe & BFS) :** Propagation de proche en proche (parcours en largeur) pour orienter les normales vers l'extérieur/l'intérieur de manière cohérente.
4. **Étalement gaussien (*Splatting*) :** Distribution de l'influence des normales sur les nœuds de la grille pour calculer le champ de divergence (second membre de l'équation de Poisson).
5. **Schéma aux différences finies :** Discrétisation de l'équation de Poisson $\Delta f = g$ par un schéma à cinq points.
6. **Solveur de Cholesky Creux (CSC) :**
   - Remplissage de la matrice au format *Compressed Sparse Column* (CSC).
   - Numérotation intelligente des nœuds (RCMK ou METIS).
   - Factorisation de Cholesky symbolique (prédiction du remplissage *fill-in*).
   - Factorisation de Cholesky numérique.
   - Résolution par substitutions avant et arrière.
7. **Extraction de la courbe :** Calcul du niveau optimal $\tau$ et interpolation pour tracer la ligne de niveau finale.

## 🛠️ Dépendances

Le code est écrit en **C**. Il requiert les bibliothèques externes suivantes :
* **FLANN** (Fast Library for Approximate Nearest Neighbors) : pour la recherche spatiale rapide des voisins.
* **METIS** *(Optionnel)* : pour calculer la permutation des nœuds visant à minimiser le remplissage lors de la factorisation.

## 🚀 Compilation et Exécution

*(À compléter avec votre Makefile ou vos commandes de compilation)*

### Exemple d'exécution

Le programme principal prend en entrée un fichier texte contenant le nuage de points (nombre de points sur la première ligne, puis coordonnées `x y`).

```bash
./reconstruction nuage.txt -N 101 -sigma 0.05 -k 10
```

**Paramètres :**
* `nuage.txt` : Chemin vers le fichier contenant le nuage de points.
* `-N` *(Obligatoire)* : Nombre de nœuds par côté de la grille (bord compris).
* `-k` *(Optionnel, défaut = 5)* : Nombre de voisins à utiliser pour l'estimation locale.
* `-sigma` *(Optionnel, défaut = 1.5 * pas de grille)* : Largeur de la gaussienne pour l'étalement.

### Sortie

Le programme génère un fichier `reconstruction.dat` (ou similaire) contenant la géométrie de la grille et les valeurs du champ scalaire $f$, prêt à être lu par un script Python ou un outil de tracé pour visualiser le résultat.

## 👥 Équipe (Groupe de 4)

* **[Asscherickx / Diego]** - *[Rôle/Tâches, ex: Solveur Cholesky]*
* **[Ledent / Mathieu]** - *[Rôle/Tâches, ex: Recherche spatiale & Normales]*
* **[Nom / Prénom Étudiant 3]** - *[Rôle/Tâches, ex: Orientation & Parcours BFS]*
* **[Nom / Prénom Étudiant 4]** - *[Rôle/Tâches, ex: Étalement gaussien & Visualisation]*


