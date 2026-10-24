"""Raffine la grille en gardant le nuage et sigma = 0.06 fixes"""
from pathlib import Path
import csv
import math
import subprocess
import sys
try:
    from .mesures_experience import ajouter_mesures
except ImportError:
    from mesures_experience import ajouter_mesures

racine = Path(__file__).resolve().parents[2]
sortie = racine / "test/resultat_experience/raffinement_grille"
sortie.mkdir(parents=True, exist_ok=True)
resume = sortie / "resume.csv"
resume.write_text("")

for taille in (41, 61, 81, 101):
    donnees = sortie / f"grille_{taille}.dat"
    image = sortie / f"grille_{taille}.png"
    subprocess.run([
        racine / "reconstruction", racine / "nuages/cercle_400.txt",
        "-N", str(taille), "-k", "10", "-sigma", "0.06",
        "-num", "nd", "-o", donnees, "-v",
    ], check=True, cwd=racine)
    subprocess.run([
        sys.executable, racine / "scripts/visualiser.py", donnees,
        "-o", image, "--reference", "cercle", "0", "0", "1",
    ], check=True, cwd=racine)
    ajouter_mesures(donnees, resume, f"grille_{taille}", 0.06, 10, cercle=True)

with resume.open() as f:
    lignes = list(csv.DictReader(f))
with (sortie / "ordres.csv").open("w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["grille_1", "grille_2", "ordre_erreur_max"])
    for a, b in zip(lignes, lignes[1:]):
        ordre = math.log(float(a["erreur_max_reference"]) / float(b["erreur_max_reference"]))
        ordre /= math.log(float(a["h"]) / float(b["h"]))
        writer.writerow([a["nx"], b["nx"], ordre])

print(f"Resultats dans {sortie}")
