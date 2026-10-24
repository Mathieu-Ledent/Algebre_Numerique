# LINMA1170 - Devoir 1 : Reconstruction de surface (Poisson 2D)

**Année académique :** 2026-2027  
**Cours :** LINMA1170 - Analyse numérique  


Ce projet implémente une chaîne complète de reconstruction de surface (courbe 2D) à partir d'un nuage de points potentiellement bruité. La méthode repose sur la **reconstruction de Poisson**, combinée à un solveur de système linéaire utilisant une **factorisation de Cholesky creuse**.

##   Objectifs du projet

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

##   Dépendances

Le code est écrit en **C**. Il requiert les bibliothèques externes suivantes :
* **FLANN** (Fast Library for Approximate Nearest Neighbors) : pour la recherche spatiale rapide des voisins.
* **METIS** *(Optionnel)* : pour calculer la permutation des nœuds visant à minimiser le remplissage lors de la factorisation.

##   Compilation et Exécution

```bash
make                 # FLANN détecté via pkg-config ; sinon recherche exhaustive O(N²) (avertissement)
make FLANN=0         # force la recherche exhaustive (référence de validation, aucune dépendance)
make METIS=1         # active la numérotation METIS_NodeND (option -num metis)
make test            # compile et lance tous les tests (build/test_*)
make clean
```

Installation des bibliothèques (Ubuntu/Debian) : `sudo apt install libflann-dev libmetis-dev`.
Visualisation : `pip install numpy matplotlib`.

### Organisation du code

| Fichier | Rôle |
|---|---|
| `src/main.c` | Programme principal (section 7.1) : lecture du nuage, grille avec marge de 10 %, export |
| `utils/adjacences.c` | k plus proches voisins avec FLANN (kd-tree) |
| `utils/adjacences_exhaustif.c` | Même fonction sans FLANN (recherche exhaustive, sélection partielle) |
| `utils/normales.c` | Droites moyennes (covariance 2×2) et orientation BFS sur le graphe **non orienté** |
| `utils/divergence.c` | Étalement gaussien : second membre `g = div n` |
| `utils/poisson.c` | Matrice `A = -Δ_h` (triangle inférieur, CSC), numérotation prise en compte à l'assemblage |
| `utils/numerotation.c` | Numérotations : ligne par ligne, RCMK, dissection emboîtée géométrique, METIS |
| `utils/cholesky.c` | Cholesky symbolique (remplissage) et numérique (left-looking) |
| `utils/solver.c` | Substitutions avant/arrière et `resoudrePoisson` |
| `utils/utils.c` | Interpolation bilinéaire, niveau τ, poids adaptés, chronomètre |
| `utils/reconstruction.c` | `reconstructionPoisson` (chaîne complète) et `exporterResultats` |
| `headers/grille.h` | Définition **unique** de la structure `Grille` |
| `scripts/generer_nuages.py` | Nuages des cas de la section 7 (cercle, concave, emboîtés, bruit, lacune...) |
| `scripts/visualiser.py` | Nuage + Γ_h, carte de f, dispersion de f_h(p_i), zoom, erreur géométrique |

`reconstructionPoisson` garde en mémoire la numérotation et le facteur L : un nouvel appel sur la
même grille (même nx, ny, h et numérotation) ne refait que l'étalement et les substitutions
(réutilisation du facteur pour plusieurs seconds membres, section 5.1).

### Exemple d'exécution

Le programme principal prend en entrée un fichier texte contenant le nuage de points (nombre de points sur la première ligne, puis coordonnées `x y`).

```bash
python3 scripts/generer_nuages.py                      # crée nuages/*.txt
./reconstruction nuages/haricot.txt -N 101 -sigma 0.05 -k 10 -num nd -v
python3 scripts/visualiser.py reconstruction.dat -o haricot.png
python3 scripts/visualiser.py reconstruction.dat --reference cercle 0 0 1 --zoom 0.6 1.1 -0.25 0.25
```

**Paramètres :**
* `nuage.txt` : Chemin vers le fichier contenant le nuage de points.
* `-N` *(Obligatoire)* : Nombre de nœuds par côté de la grille (bord compris).
* `-k` *(Optionnel, défaut = 5)* : Nombre de voisins à utiliser pour l'estimation locale.
* `-sigma` *(Optionnel, défaut = 1.5 * pas de grille)* : Largeur de la gaussienne pour l'étalement.
* `-num ligne|rcmk|nd|metis` *(Optionnel, défaut = ligne)* : Numérotation des nœuds intérieurs.
* `-poids` *(Optionnel)* : poids a_i proportionnels à l'espacement local (échantillonnage non uniforme).
* `-inverser tout|x` *(Optionnel)* : inverse toutes les normales ou seulement leur composante x (étude 6.2).
* `-o fichier` *(Optionnel, défaut = reconstruction.dat)* : fichier de sortie.
* `-v` *(Optionnel)* : affiche les temps de chaque étape, nnz(A), nnz(L), le résidu et la dispersion autour de τ.

### Sortie

Le programme génère `reconstruction.dat` : le nuage (`N` puis `x y`), la grille (`nx ny xmin ymin h`),
le niveau τ, puis f en `ny` lignes de `nx` valeurs (ligne s : `f[r + nx*s]`). Les blocs sont
séparés par des lignes de commentaire `#`. `scripts/visualiser.py` lit ce format.

### Quelques résultats (`make test`)

* Solution fabriquée `sin(πx) sin(πy)` : erreur divisée par 4 quand h est divisé par 2 (ordre 2.00), résidu relatif ~1e-12.
* Remplissage pour nx = 257 (M = 65 025) : nnz(L) = 16.6 M (ligne), 11.2 M (RCMK), 2.3 M (dissection emboîtée) ;
  factorisation numérique 1.0 s / 0.50 s / 0.07 s. Un second membre supplémentaire ne coûte que les substitutions.

##   Équipe (Groupe de 4)

* **[Asscherickx / Diego]** - *[Résoluton par Cholesky creux et visualisation]*
* **[Ledent /Mathieu]** - *[Orientation et parcours BFS + Schema de poisson]*
* **[Vos de Wael / Thibault]** - *[Recherche spatiale et estimation des normales]*
* **[Herygers / Maximilien]** - *[Étalement gaussien et expériences]*

