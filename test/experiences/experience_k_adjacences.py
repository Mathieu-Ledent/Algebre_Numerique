"""Fait varier k sur le haricot et montre le graphe des voisins"""
from pathlib import Path
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np

try:
    from .mesures_experience import ajouter_mesures, extraire_residu
except ImportError:
    from mesures_experience import ajouter_mesures, extraire_residu

racine = Path(__file__).resolve().parents[2]
sortie = racine / "test/resultat_experience/k_adjacences"
sortie.mkdir(parents=True, exist_ok=True)
resume = sortie / "resume.csv"
resume.write_text("")

nuage = racine / "nuages/haricot.txt"
points = np.loadtxt(nuage, skiprows=1)
d2 = ((points[:, None, :] - points[None, :, :]) ** 2).sum(axis=2)
valeurs_k = (5, 10, 20, 40)
fig, axes = plt.subplots(2, 2, figsize=(9, 9))

for k, ax in zip(valeurs_k, axes.flat):
    voisins = np.argsort(d2, axis=1)[:, 1:k + 1]
    aretes = np.stack((np.repeat(points, k, axis=0), points[voisins.ravel()]), axis=1)
    ax.add_collection(LineCollection(aretes, colors="C0", linewidths=0.25, alpha=0.08))
    ax.plot(points[:, 0], points[:, 1], ".", color="black", ms=1)
    ax.set_title(f"k = {k}")
    ax.set_aspect("equal")
    ax.autoscale()

    donnees = sortie / f"k_{k}.dat"
    image = sortie / f"k_{k}.png"
    calcul = subprocess.run([
        racine / "reconstruction", nuage, "-N", "81", "-k", str(k),
        "-sigma", "0.06", "-num", "nd", "-o", donnees, "-v",
    ], check=True, cwd=racine, capture_output=True, text=True)
    print(calcul.stdout, end="")
    subprocess.run([
        sys.executable, racine / "scripts/visualiser.py", donnees,
        "-o", image, "--zoom", "-0.7", "0.7", "-0.9", "0.9",
    ], check=True, cwd=racine)
    ajouter_mesures(donnees, resume, f"k_{k}", 0.06, k,
                    residu=extraire_residu(calcul.stdout))

fig.tight_layout()
fig.savefig(sortie / "aretes_knn.png", dpi=150)
plt.close(fig)
print(f"Resultats dans {sortie}")
