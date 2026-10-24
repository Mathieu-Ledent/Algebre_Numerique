"""Mesures simples pour les scripts d'experience"""
import csv
from pathlib import Path
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

racine = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(racine / "scripts"))
from visualiser import bilineaire, distances_points_segments, lire


def extraire_residu(texte):
    resultat = re.search(r"Residu relatif .* = ([0-9.eE+-]+)", texte)
    return float(resultat.group(1)) if resultat else np.nan


def ajouter_mesures(donnees, fichier_csv, cas, sigma, k, poids=False,
                    cercle=False, residu=np.nan):
    points, nx, ny, xmin, ymin, h, tau, f = lire(donnees)
    x = xmin + h * np.arange(nx)
    y = ymin + h * np.arange(ny)
    fp = bilineaire(f, xmin, ymin, h, points)

    fig, ax = plt.subplots()
    contour = ax.contour(x, y, f, levels=[tau])
    segments = [s for s in contour.allsegs[0] if len(s) > 1]
    plt.close(fig)

    longueur = sum(np.sqrt(((s[1:] - s[:-1]) ** 2).sum(axis=1)).sum() for s in segments)
    aire = 0.0
    for s in segments:
        aire += 0.5 * abs(np.sum(s[:, 0] * np.roll(s[:, 1], -1)
                                - s[:, 1] * np.roll(s[:, 0], -1)))
    if segments:
        erreur_points = distances_points_segments(points, segments)
        sommets = np.vstack(segments)
        erreur_sommets = [np.sqrt(((points - sommet) ** 2).sum(axis=1)).min()
                           for sommet in sommets]
    else:
        erreur_points = np.array([np.inf])
        erreur_sommets = [np.inf]

    voisins = []
    for i, point in enumerate(points):
        d2 = ((points - point) ** 2).sum(axis=1)
        d2[i] = np.inf
        voisins.append(np.sqrt(d2.min()))
    moyen = float(np.mean(voisins))
    median = float(np.median(voisins))

    emax1 = emoy1 = emax2 = emoy2 = np.nan
    if cercle:
        echantillons = np.vstack(segments)
        e1 = np.abs(np.hypot(echantillons[:, 0], echantillons[:, 1]) - 1.0)
        t = np.linspace(0, 2 * np.pi, 4000)
        reference = np.c_[np.cos(t), np.sin(t)]
        e2 = distances_points_segments(reference, segments)
        emax1, emoy1 = float(e1.max()), float(e1.mean())
        emax2, emoy2 = float(e2.max()), float(e2.mean())

    colonnes = [
        "cas", "points", "nx", "h", "k", "sigma", "sigma_sur_h",
        "espacement_moyen", "espacement_median", "sigma_sur_espacement",
        "poids", "tau", "ecart_type_niveau", "composantes", "fermees",
        "longueur_contour", "aire_contour", "erreur_max_nuage_contour",
        "erreur_moy_nuage_contour", "erreur_max_contour_nuage",
        "erreur_moy_contour_nuage",
        "erreur_max_reconstruction", "erreur_moy_reconstruction",
        "erreur_max_reference", "erreur_moy_reference", "residu_relatif",
    ]
    ligne = [
        cas, len(points), nx, h, k, sigma, sigma / h, moyen, median,
        sigma / median, int(poids), tau, float(np.std(fp - tau)), len(segments),
        sum(np.allclose(s[0], s[-1]) for s in segments),
        longueur, aire, float(erreur_points.max()), float(erreur_points.mean()),
        float(np.max(erreur_sommets)), float(np.mean(erreur_sommets)),
        emax1, emoy1, emax2, emoy2, residu,
    ]

    fichier_csv = Path(fichier_csv)
    nouveau = not fichier_csv.exists() or fichier_csv.stat().st_size == 0
    with fichier_csv.open("a", newline="") as fcsv:
        writer = csv.writer(fcsv)
        if nouveau:
            writer.writerow(colonnes)
        writer.writerow(ligne)
