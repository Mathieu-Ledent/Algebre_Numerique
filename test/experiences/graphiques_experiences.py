"""Cree les graphes de synthese a partir des CSV des experiences"""
from pathlib import Path
import csv

import matplotlib.pyplot as plt
import numpy as np

racine = Path(__file__).resolve().parents[1]
base = racine / "resultat_experience"
sortie = base / "graphiques"
sortie.mkdir(parents=True, exist_ok=True)


def lire(nom):
    with (base / nom / "resume.csv").open() as f:
        return list(csv.DictReader(f))


def valeurs(lignes, nom):
    return np.array([float(ligne[nom]) for ligne in lignes])


normalisation = list(csv.DictReader((base / "normalisation_gaussienne/normalisation.csv").open()))
x = valeurs(normalisation, "sigma_sur_h")
plt.figure()
plt.semilogy(x, valeurs(normalisation, "erreur"), "o-")
plt.xlabel(r"$\sigma/h$")
plt.ylabel(r"$|h^2\sum G_\sigma-1|$")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(sortie / "normalisation.png", dpi=150)
plt.close()

bruit = lire("cercle_bruite_sigma")
x = valeurs(bruit, "sigma_sur_h")
plt.figure()
plt.plot(x, valeurs(bruit, "erreur_max_reference"), "o-", label="erreur maximale")
plt.plot(x, valeurs(bruit, "erreur_moy_reference"), "o-", label="erreur moyenne")
plt.xlabel(r"$\sigma/h$")
plt.ylabel("erreur geometrique")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(sortie / "cercle_bruite_erreur.png", dpi=150)
plt.close()

plt.figure()
plt.plot(x, valeurs(bruit, "ecart_type_niveau"), "o-")
plt.xlabel(r"$\sigma/h$")
plt.ylabel(r"ecart-type de $f_h(p_i)-\tau$")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(sortie / "cercle_bruite_dispersion.png", dpi=150)
plt.close()

poids = lire("non_uniforme_poids")
labels = ["sans poids", "avec poids"]
positions = np.arange(len(labels))
fig, axes = plt.subplots(1, 2, figsize=(9, 4))
axes[0].bar(positions - 0.18, valeurs(poids, "erreur_max_reference"), 0.36, label="max")
axes[0].bar(positions + 0.18, valeurs(poids, "erreur_moy_reference"), 0.36, label="moyenne")
axes[0].set_xticks(positions, labels)
axes[0].set_ylabel("erreur geometrique")
axes[0].legend()
axes[1].bar(labels, valeurs(poids, "ecart_type_niveau"))
axes[1].set_ylabel(r"ecart-type de $f_h(p_i)-\tau$")
fig.tight_layout()
fig.savefig(sortie / "poids_non_uniformes.png", dpi=150)
plt.close(fig)

grilles = lire("raffinement_grille")
h = valeurs(grilles, "h")
emax = valeurs(grilles, "erreur_max_reference")
emoy = valeurs(grilles, "erreur_moy_reference")
plt.figure()
plt.loglog(h, emax, "o-", label="erreur maximale")
plt.loglog(h, emoy, "o-", label="erreur moyenne")
plt.loglog(h, emax[0] * (h / h[0]) ** 2, "--", label=r"reference $h^2$")
plt.xlabel("pas h")
plt.ylabel("erreur geometrique")
plt.grid(alpha=0.3, which="both")
plt.legend()
plt.tight_layout()
plt.savefig(sortie / "raffinement_grille.png", dpi=150)
plt.close()

nuages = lire("raffinement_nuage")
n = valeurs(nuages, "points")
plt.figure()
plt.semilogx(n, valeurs(nuages, "erreur_max_reference"), "o-", label="erreur maximale")
plt.semilogx(n, valeurs(nuages, "erreur_moy_reference"), "o-", label="erreur moyenne")
plt.xlabel("nombre de points")
plt.ylabel("erreur geometrique")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(sortie / "raffinement_nuage.png", dpi=150)
plt.close()

haricot = lire("haricot_sigma")
x = valeurs(haricot, "sigma_sur_h")
erreur = np.maximum(valeurs(haricot, "erreur_max_nuage_contour"),
                    valeurs(haricot, "erreur_max_contour_nuage"))
fig, ax1 = plt.subplots()
ax1.plot(x, valeurs(haricot, "longueur_contour"), "o-", color="C0")
ax1.set_xlabel(r"$\sigma/h$")
ax1.set_ylabel("longueur du contour", color="C0")
ax2 = ax1.twinx()
ax2.plot(x, erreur, "s-", color="C3")
ax2.set_ylabel("distance maximale au nuage", color="C3")
ax1.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(sortie / "haricot_lissage.png", dpi=150)
plt.close(fig)

print(f"Graphes dans {sortie}")
