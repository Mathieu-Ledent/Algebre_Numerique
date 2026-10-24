"""Verifie que h^2 somme G_sigma approche 1 quand la grille resout la gaussienne"""
from pathlib import Path
import csv

import matplotlib.pyplot as plt
import numpy as np

racine = Path(__file__).resolve().parents[1]
sortie = racine / "resultat_experience" / "normalisation_gaussienne"
sortie.mkdir(parents=True, exist_ok=True)

sigma = 1.0
rapports = (0.5, 1.0, 1.5, 2.0, 4.0)
sommes = []

for rapport in rapports:
    h = sigma / rapport
    x = np.arange(-6 * sigma, 6 * sigma + h / 2, h)
    X, Y = np.meshgrid(x, x)
    G = np.exp(-(X * X + Y * Y) / (2 * sigma * sigma))
    G /= 2 * np.pi * sigma * sigma
    sommes.append(h * h * G.sum())

with (sortie / "normalisation.csv").open("w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["sigma_sur_h", "somme", "erreur"])
    for rapport, somme in zip(rapports, sommes):
        writer.writerow([rapport, somme, abs(somme - 1)])

plt.plot(rapports, np.abs(np.array(sommes) - 1), "o-")
plt.yscale("log")
plt.xlabel(r"$\sigma/h$")
plt.ylabel(r"$|h^2\sum G_\sigma-1|$")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(sortie / "normalisation.png", dpi=150)
print(f"Resultats dans {sortie}")
