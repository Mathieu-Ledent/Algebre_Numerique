"""Teste deux composantes et un anneau avec trois ordres de points"""
from pathlib import Path
import subprocess
import sys

import numpy as np

try:
    from .mesures_experience import ajouter_mesures, extraire_residu
except ImportError:
    from mesures_experience import ajouter_mesures, extraire_residu

racine = Path(__file__).resolve().parents[2]
sortie = racine / "test/resultat_experience/composantes"
sortie.mkdir(parents=True, exist_ok=True)
resume = sortie / "resume.csv"
resume.write_text("")

for cas in ("deux_cercles", "anneau"):
    points = np.loadtxt(racine / f"nuages/{cas}.txt", skiprows=1)
    for ordre in (1, 2, 3):
        melanges = np.random.default_rng(ordre).permutation(points)
        nuage = sortie / f"{cas}_ordre_{ordre}.txt"
        np.savetxt(nuage, melanges, header=str(len(melanges)), comments="")
        nom = f"{cas}_ordre_{ordre}"
        donnees = sortie / f"{nom}.dat"
        image = sortie / f"{nom}.png"

        calcul = subprocess.run([
            racine / "reconstruction", nuage, "-N", "101", "-k", "10",
            "-sigma", "0.06", "-num", "nd", "-o", donnees, "-v",
        ], check=True, cwd=racine, capture_output=True, text=True)
        print(calcul.stdout, end="")
        subprocess.run([
            sys.executable, racine / "scripts/visualiser.py", donnees, "-o", image,
        ], check=True, cwd=racine)
        ajouter_mesures(donnees, resume, nom, 0.06, 10,
                        residu=extraire_residu(calcul.stdout))

print(f"Resultats dans {sortie}")
