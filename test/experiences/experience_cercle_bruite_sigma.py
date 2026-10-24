"""Reconstruit un cercle bruite avec cinq largeurs gaussiennes"""
from pathlib import Path
import subprocess
import sys
try:
    from .mesures_experience import ajouter_mesures
except ImportError:
    from mesures_experience import ajouter_mesures

racine = Path(__file__).resolve().parents[2]
sortie = racine / "test/resultat_experience/cercle_bruite_sigma"
sortie.mkdir(parents=True, exist_ok=True)
resume = sortie / "resume.csv"
resume.write_text("")

for sigma in (0.02, 0.04, 0.08, 0.14, 0.24):
    nom = f"sigma_{sigma:.2f}".replace(".", "p")
    donnees = sortie / f"{nom}.dat"
    image = sortie / f"{nom}.png"
    subprocess.run([
        racine / "reconstruction", racine / "nuages/cercle_bruit_0.06.txt",
        "-N", "81", "-k", "10", "-sigma", str(sigma),
        "-num", "nd", "-o", donnees, "-v",
    ], check=True, cwd=racine)
    subprocess.run([
        sys.executable, racine / "scripts/visualiser.py", donnees,
        "-o", image, "--reference", "cercle", "0", "0", "1",
    ], check=True, cwd=racine)
    ajouter_mesures(donnees, resume, nom, sigma, 10, cercle=True)

print(f"Resultats dans {sortie}")
