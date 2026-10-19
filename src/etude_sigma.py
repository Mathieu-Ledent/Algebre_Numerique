import numpy as np
import matplotlib.pyplot as plt


def creer_carre(n_par_cote=40):

    t = np.linspace(-0.5, 0.5, n_par_cote, endpoint=False)

    points = np.vstack((
        np.column_stack((t, -0.5 * np.ones_like(t))),
        np.column_stack((0.5 * np.ones_like(t), t)),
        np.column_stack((-t, 0.5 * np.ones_like(t))),
        np.column_stack((-0.5 * np.ones_like(t), -t)),
    ))
    normales = np.repeat(
        [[0, -1], [1, 0], [0, 1], [-1, 0]], n_par_cote, axis=0
    )
    return points, normales


def calculer_divergence(points, normales, sigma, X, Y):

    divergence = np.zeros_like(X)

    for point, normale in zip(points, normales):
        dx = X - point[0]
        dy = Y - point[1]
        distance2 = dx**2 + dy**2

        gaussienne = np.exp(-distance2 / (2 * sigma**2))
        gaussienne /= 2 * np.pi * sigma**2
        produit_scalaire = normale[0] * dx + normale[1] * dy

        divergence += produit_scalaire * gaussienne

    return -divergence / sigma**2

x = np.linspace(-0.7, 0.7, 121)
y = np.linspace(-0.7, 0.7, 121)
X, Y = np.meshgrid(x, y)
h = x[1] - x[0]

points, normales = creer_carre()
generateur = np.random.default_rng(4)
direction_bruit = generateur.normal(size=points.shape)

sigmas = [1 * h, 1.5 * h, 2 * h, 2.5 * h, 3 * h]
bruits = [0, 0.5 * h, 2 * h]

fig, axes = plt.subplots(len(bruits), len(sigmas), figsize=(10, 9))

for ligne, bruit in enumerate(bruits):
    points_bruites = points + bruit * direction_bruit

    for colonne, sigma in enumerate(sigmas):
        div = calculer_divergence(points_bruites, normales, sigma, X, Y)
        div /= np.max(np.abs(div))  

        ax = axes[ligne, colonne]
        image = ax.imshow(
            div,
            extent=[x.min(), x.max(), y.min(), y.max()],
            origin="lower",
            cmap="RdBu_r",
            vmin=-1,
            vmax=1,
        )
        ax.scatter(points_bruites[:, 0], points_bruites[:, 1], s=2, c="black")
        ax.set_title(rf"$\sigma = {sigma / h:.2g}h$")
        ax.set_ylabel(rf"bruit = {bruit / h:.2g}h")
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])

fig.colorbar(image, ax=axes, shrink=0.75, label="divergence normalisee")
fig.suptitle("Influence de sigma sur la divergence d'un carre bruite")
plt.savefig("comparaison_sigma.png", dpi=200, bbox_inches="tight")
plt.show()
