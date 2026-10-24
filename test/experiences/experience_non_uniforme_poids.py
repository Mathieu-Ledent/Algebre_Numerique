"""Compare la reconstruction non uniforme avec et sans poids adaptes"""
from pathlib import Path
import subprocess
import sys
try:
    from .mesures_experience import ajouter_mesures
except ImportError:
    from mesures_experience import ajouter_mesures

racine = Path(__file__).resolve().parents[2]
sortie = racine / "test/resultat_experience/non_uniforme_poids"
sortie.mkdir(parents=True, exist_ok=True)
resume = sortie / "resume.csv"
resume.write_text("")

for nom, poids_adaptes in (("sans_poids", False), ("avec_poids", True)):
    donnees = sortie / f"{nom}.dat"
    image = sortie / f"{nom}.png"
    commande = [
        racine / "reconstruction", racine / "nuages/non_uniforme.txt",
        "-N", "81", "-k", "10", "-sigma", "0.06",
        "-num", "nd", "-o", donnees, "-v",
    ]
    if not poids_adaptes:
        commande.append("-poids-unitaires")
    subprocess.run(commande, check=True, cwd=racine)
    subprocess.run([
        sys.executable, racine / "scripts/visualiser.py", donnees,
        "-o", image, "--reference", "cercle", "0", "0", "1",
    ], check=True, cwd=racine)
    ajouter_mesures(donnees, resume, nom, 0.06, 10, poids_adaptes, cercle=True)

print(f"Resultats dans {sortie}")
