"""Lance toutes les experiences puis cree les graphes de synthese"""
from pathlib import Path
import subprocess
import sys

dossier = Path(__file__).resolve().parent

scripts = (
    "experience_haricot_sigma.py",
    "experience_cercle_bruite_sigma.py",
    "experience_non_uniforme_poids.py",
    "experience_raffinement_grille.py",
    "experience_raffinement_nuage.py",
    "experience_normalisation_gaussienne.py",
    "experience_k_adjacences.py",
    "experience_composantes.py",
    "experience_cas_limites.py",
    "graphiques_experiences.py",
)

for script in scripts:
    print(f"===== {script}", flush=True)
    subprocess.run([sys.executable, dossier / script], check=True)
