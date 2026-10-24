"""Raffine le nuage en gardant la grille et sigma = 0.06 fixes"""
from pathlib import Path
import subprocess
import sys
try:
    from .mesures_experience import ajouter_mesures
except ImportError:
    from mesures_experience import ajouter_mesures

racine = Path(__file__).resolve().parents[2]
sortie = racine / "test/resultat_experience/raffinement_nuage"
sortie.mkdir(parents=True, exist_ok=True)
resume = sortie / "resume.csv"
resume.write_text("")

for points in (100, 400, 1600):
    nom = f"cercle_{points}"
    donnees = sortie / f"{nom}.dat"
    image = sortie / f"{nom}.png"
    subprocess.run([
        racine / "reconstruction", racine / f"nuages/{nom}.txt",
        "-N", "101", "-k", "10", "-sigma", "0.06",
        "-num", "nd", "-o", donnees, "-v",
    ], check=True, cwd=racine)
    subprocess.run([
        sys.executable, racine / "scripts/visualiser.py", donnees,
        "-o", image, "--reference", "cercle", "0", "0", "1",
    ], check=True, cwd=racine)
    ajouter_mesures(donnees, resume, nom, 0.06, 10, cercle=True)

print(f"Resultats dans {sortie}")
