#!/usr/bin/env python3
"""Genere les nuages de points des cas etudies (section 7 de l'enonce) dans le dossier nuages/.

Format : nombre de points sur la premiere ligne, puis une ligne "x y" par point.
Bibliotheque standard uniquement ; graine fixe pour des configurations reproductibles.

    python3 scripts/generer_nuages.py [dossier]
"""
import math
import os
import random
import sys


def cercle(n, cx=0.0, cy=0.0, r=1.0, bruit=0.0, sens=1):
    pts = []
    for i in range(n):
        t = sens * 2 * math.pi * i / n
        rr = r + random.gauss(0.0, bruit) if bruit > 0 else r
        pts.append((cx + rr * math.cos(t), cy + rr * math.sin(t)))
    return pts


def ellipse(n, a=1.5, b=0.7, bruit=0.0):
    return [(a * math.cos(t) + random.gauss(0, bruit), b * math.sin(t) + random.gauss(0, bruit))
            for t in (2 * math.pi * i / n for i in range(n))]


def courbe_polaire(n, rayon):
    """Courbe fermee r = rayon(t), echantillonnee uniformement en t."""
    return [(rayon(t) * math.cos(t), rayon(t) * math.sin(t))
            for t in (2 * math.pi * i / n for i in range(n))]


def haricot(n):
    """Courbe fortement concave : deux lobes relies par un col etroit (branches proches)."""
    return courbe_polaire(n, lambda t: 1.0 + 0.75 * math.cos(2 * t) ** 3)


def etoile(n, branches=5, creux=0.55):
    return courbe_polaire(n, lambda t: 1.0 - creux * 0.5 * (1 - math.cos(branches * t)))


def non_uniforme(n, bruit=0.02):
    """Cercle echantillonne de facon tres non uniforme (densite x10 d'un cote) et bruite."""
    pts = []
    while len(pts) < n:
        t = random.uniform(0, 2 * math.pi)
        if random.random() < 0.1 + 0.9 * (0.5 + 0.5 * math.cos(t)) ** 4:
            r = 1.0 + random.gauss(0, bruit)
            pts.append((r * math.cos(t), r * math.sin(t)))
    return pts


def arc_ouvert(n, ouverture=math.pi / 2):
    """Arc de cercle ouvert (lacune d'angle 'ouverture')."""
    debut, fin = ouverture / 2, 2 * math.pi - ouverture / 2
    return [(math.cos(t), math.sin(t)) for t in
            (debut + (fin - debut) * i / (n - 1) for i in range(n))]


def avec_aberrants(pts, m, etendue=1.3):
    return pts + [(random.uniform(-etendue, etendue), random.uniform(-etendue, etendue)) for _ in range(m)]


def ecrire(dossier, nom, pts):
    chemin = os.path.join(dossier, nom + ".txt")
    random.shuffle(pts)            # l'ordre des points dans le fichier ne doit pas compter
    with open(chemin, "w") as f:
        f.write(f"{len(pts)}\n")
        for x, y in pts:
            f.write(f"{x:.12g} {y:.12g}\n")
    print(f"{chemin:40s} {len(pts):6d} points")


def main():
    dossier = sys.argv[1] if len(sys.argv) > 1 else "nuages"
    os.makedirs(dossier, exist_ok=True)
    random.seed(1170)

    # Validation de reference et convergence (nuage et grille raffines).
    for n in (100, 400, 1600):
        ecrire(dossier, f"cercle_{n}", cercle(n))
    ecrire(dossier, "ellipse", ellipse(500))
    # Courbes concaves et branches proches.
    ecrire(dossier, "haricot", haricot(800))
    ecrire(dossier, "etoile", etoile(1000))
    # Deux composantes separees ; deux contours emboites (anneau, domaine avec trou).
    ecrire(dossier, "deux_cercles", cercle(300, -1.2, 0, 0.8) + cercle(300, 1.2, 0, 0.8))
    ecrire(dossier, "anneau", cercle(400, r=1.0) + cercle(250, r=0.5))
    # Bruit et densite non uniforme.
    for b in (0.01, 0.03, 0.06):
        ecrire(dossier, f"cercle_bruit_{b:g}", cercle(600, bruit=b))
    ecrire(dossier, "non_uniforme", non_uniforme(600))
    # Lacune importante, arc ouvert, points aberrants.
    ecrire(dossier, "lacune", arc_ouvert(400, math.pi / 3))
    ecrire(dossier, "arc_ouvert", arc_ouvert(300, math.pi))
    ecrire(dossier, "aberrants", avec_aberrants(cercle(500, bruit=0.01), 25))


if __name__ == "__main__":
    main()
