#!/usr/bin/env python3
"""Trace le resultat de ./reconstruction (fichier reconstruction.dat).

Figure : (1) nuage + ligne de niveau Gamma_h (rapport d'aspect egal), (2) carte de f avec le
niveau tau, (3) dispersion des f_h(p_i) autour de tau. La ligne de niveau est extraite par
matplotlib.contour (algorithme "marching squares" de la bibliotheque contourpy : interpolation
lineaire sur les aretes des cellules ; les cellules ambigues (point selle) sont tranchees par
la valeur moyenne des quatre coins, comportement par defaut de contourpy).

    python3 scripts/visualiser.py reconstruction.dat [-o figure.png] [--zoom x0 x1 y0 y1]
                                  [--reference cercle CX CY R] [--niveau TAU]

Dependances : numpy, matplotlib.
"""
import argparse
import sys

import numpy as np


def lire(chemin):
    with open(chemin) as f:
        lignes = [l for l in f if not l.startswith("#")]
    n = int(lignes[0])
    points = np.loadtxt(lignes[1:1 + n], ndmin=2)
    nx, ny, xmin, ymin, h = lignes[1 + n].split()
    nx, ny, xmin, ymin, h = int(nx), int(ny), float(xmin), float(ymin), float(h)
    tau = float(lignes[2 + n])
    f = np.loadtxt(lignes[3 + n:3 + n + ny], ndmin=2)   # f[s, r] = f_{r,s}
    assert f.shape == (ny, nx), "dimensions de f incoherentes"
    return points, nx, ny, xmin, ymin, h, tau, f


def bilineaire(f, xmin, ymin, h, pts):
    ny, nx = f.shape
    u = np.clip((pts[:, 0] - xmin) / h, 0, nx - 1)
    v = np.clip((pts[:, 1] - ymin) / h, 0, ny - 1)
    r = np.minimum(np.floor(u).astype(int), nx - 2)
    s = np.minimum(np.floor(v).astype(int), ny - 2)
    a, b = u - r, v - s
    return ((1 - a) * (1 - b) * f[s, r] + a * (1 - b) * f[s, r + 1]
            + (1 - a) * b * f[s + 1, r] + a * b * f[s + 1, r + 1])


def distances_points_segments(p, segments):
    """Distance de chaque point de p (m x 2) a la polyligne la plus proche (liste de k x 2)."""
    d = np.full(len(p), np.inf)
    for seg in segments:
        if len(seg) < 2:
            continue
        a, b = seg[:-1], seg[1:]
        ab = b - a
        l2 = np.maximum((ab ** 2).sum(1), 1e-300)
        for i0 in range(0, len(p), 2000):           # par blocs pour limiter la memoire
            q = p[i0:i0 + 2000, None, :]
            t = np.clip(((q - a) * ab).sum(2) / l2, 0, 1)
            proj = a + t[..., None] * ab
            d[i0:i0 + 2000] = np.minimum(d[i0:i0 + 2000], np.sqrt(((q - proj) ** 2).sum(2)).min(1))
    return d


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fichier", nargs="?", default="reconstruction.dat")
    ap.add_argument("-o", "--sortie", help="enregistre la figure au lieu de l'afficher")
    ap.add_argument("--zoom", nargs=4, type=float, metavar=("X0", "X1", "Y0", "Y1"),
                    help="ajoute un agrandissement sur une zone difficile")
    ap.add_argument("--reference", nargs=4, metavar=("cercle", "CX", "CY", "R"),
                    help="courbe exacte (cas synthetique) pour l'erreur geometrique")
    ap.add_argument("--niveau", type=float, help="trace un autre niveau que tau (sensibilite)")
    args = ap.parse_args()

    import matplotlib
    if args.sortie:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    points, nx, ny, xmin, ymin, h, tau, f = lire(args.fichier)
    niveau = tau if args.niveau is None else args.niveau
    x = xmin + h * np.arange(nx)
    y = ymin + h * np.arange(ny)

    fp = bilineaire(f, xmin, ymin, h, points)
    print(f"Grille {nx} x {ny}, h = {h:.4g}, {len(points)} points")
    print(f"tau = {tau:.6g} ; f sur la grille dans [{f.min():.4g}, {f.max():.4g}]")
    print(f"f_h(p_i) - tau : ecart-type {np.std(fp - tau):.3g}, "
          f"min {np.min(fp - tau):.3g}, max {np.max(fp - tau):.3g}")
    if abs(f.max() - f.min()) > 0:
        print(f"ecart-type relatif a l'amplitude de f : {np.std(fp - tau) / (f.max() - f.min()):.3g}")

    nb = 3 + (args.zoom is not None)
    fig, axes = plt.subplots(1, nb, figsize=(5.2 * nb, 5))

    ax = axes[0]
    cs = ax.contour(x, y, f, levels=[niveau], colors="C3", linewidths=1.5)
    segments = [s for s in cs.allsegs[0] if len(s) > 1]
    fermees = sum(np.allclose(s[0], s[-1]) for s in segments)
    print(f"Gamma_h : {len(segments)} composante(s), dont {fermees} fermee(s) ; "
          f"longueurs (nb sommets) : {sorted((len(s) for s in segments), reverse=True)}")
    ax.plot(points[:, 0], points[:, 1], ".", ms=2.5, color="0.2", label="nuage")
    ax.plot([], [], color="C3", label=f"$f_h = {niveau:.4g}$")
    ax.set_aspect("equal")
    ax.set_title("Nuage et courbe reconstruite")

    if args.reference:
        cx, cy, r = map(float, args.reference[1:])
        t = np.linspace(0, 2 * np.pi, 4000)
        ref = np.c_[cx + r * np.cos(t), cy + r * np.sin(t)]
        ax.plot(ref[:, 0], ref[:, 1], "--", color="C2", lw=0.8, label="reference")
        ech = np.vstack(segments) if segments else np.zeros((0, 2))
        d1 = np.abs(np.hypot(ech[:, 0] - cx, ech[:, 1] - cy) - r) if len(ech) else np.array([np.inf])
        d2 = distances_points_segments(ref, segments) if segments else np.array([np.inf])
        print(f"Erreur geometrique : reconstruite -> reference max {d1.max():.4g} (moy {d1.mean():.4g}) ; "
              f"reference -> reconstruite max {d2.max():.4g} (moy {d2.mean():.4g}) ; h = {h:.4g}")
    axes[0].legend(loc="upper right", fontsize=8)

    ax = axes[1]
    im = ax.pcolormesh(x, y, f, shading="nearest", cmap="RdBu_r")
    ax.contour(x, y, f, levels=[niveau], colors="k", linewidths=1)
    fig.colorbar(im, ax=ax, shrink=0.8, label="f")
    ax.set_aspect("equal")
    ax.set_title(f"Carte de f, niveau $\\tau$ = {niveau:.4g}")

    ax = axes[2]
    ax.hist(fp - tau, bins=40, color="C0")
    ax.axvline(0, color="C3")
    ax.set_xlabel(r"$f_h(p_i) - \tau$")
    ax.set_title("Dispersion autour du niveau")

    if args.zoom:
        x0, x1, y0, y1 = args.zoom
        ax = axes[3]
        ax.contour(x, y, f, levels=[niveau], colors="C3", linewidths=1.5)
        ax.plot(points[:, 0], points[:, 1], ".", ms=4, color="0.2")
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        ax.set_aspect("equal")
        ax.set_xticks(x[(x >= x0) & (x <= x1)], minor=True)
        ax.set_yticks(y[(y >= y0) & (y <= y1)], minor=True)
        ax.grid(which="minor", lw=0.3)
        ax.set_title("Agrandissement (grille en fond)")

    fig.tight_layout()
    if args.sortie:
        fig.savefig(args.sortie, dpi=150)
        print(f"Figure enregistree : {args.sortie}")
    else:
        plt.show()


if __name__ == "__main__":
    sys.exit(main())