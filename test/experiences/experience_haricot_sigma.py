"""Reconstruit le haricot avec cinq largeurs gaussiennes"""
from pathlib import Path
import subprocess
import sys
try:
    from .mesures_experience import ajouter_mesures
except ImportError:
    from mesures_experience import ajouter_mesures

racine = Path(__file__).resolve().parents[2]
sortie = racine / "test/resultat_experience/haricot_sigma"
sortie.mkdir(parents=True, exist_ok=True)
resume = sortie / "resume.csv"
resume.write_text("")

for sigma in (0.03, 0.06, 0.12, 0.24, 0.40):
    nom = f"sigma_{sigma:.2f}".replace(".", "p")
    donnees = sortie / f"{nom}.dat"
    image = sortie / f"{nom}.png"
    subprocess.run([
        racine / "reconstruction", racine / "nuages/haricot.txt",
        "-N", "81", "-k", "10", "-sigma", str(sigma),
        "-num", "nd", "-o", donnees, "-v",
    ], check=True, cwd=racine)
    subprocess.run([
        sys.executable, racine / "scripts/visualiser.py", donnees,
        "-o", image, "--zoom", "-0.7", "0.7", "-0.9", "0.9",
    ], check=True, cwd=racine)
    ajouter_mesures(donnees, resume, nom, sigma, 10)

print(f"Resultats dans {sortie}")
