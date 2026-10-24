"""Reconstruit les nuages avec lacune, arc ouvert et points aberrants"""
from pathlib import Path
import subprocess
import sys

try:
    from .mesures_experience import ajouter_mesures, extraire_residu
except ImportError:
    from mesures_experience import ajouter_mesures, extraire_residu

racine = Path(__file__).resolve().parents[2]
sortie = racine / "test/resultat_experience/cas_limites"
sortie.mkdir(parents=True, exist_ok=True)
resume = sortie / "resume.csv"
resume.write_text("")

for cas in ("lacune", "arc_ouvert", "aberrants"):
    donnees = sortie / f"{cas}.dat"
    image = sortie / f"{cas}.png"
    calcul = subprocess.run([
        racine / "reconstruction", racine / f"nuages/{cas}.txt",
        "-N", "81", "-k", "10", "-sigma", "0.06",
        "-num", "nd", "-o", donnees, "-v",
    ], check=True, cwd=racine, capture_output=True, text=True)
    print(calcul.stdout, end="")
    subprocess.run([
        sys.executable, racine / "scripts/visualiser.py", donnees, "-o", image,
    ], check=True, cwd=racine)
    ajouter_mesures(donnees, resume, cas, 0.06, 10,
                    residu=extraire_residu(calcul.stdout))

print(f"Resultats dans {sortie}")
