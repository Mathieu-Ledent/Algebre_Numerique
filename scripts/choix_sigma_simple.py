import numpy as np
import matplotlib.pyplot as plt


def creer_carre(n=40):

    t = np.linspace(-0.5, 0.5, n, endpoint=False)

    points = np.vstack([
        np.column_stack([t, -0.5 * np.ones(n)]),
        np.column_stack([0.5 * np.ones(n), t]),
        np.column_stack([-t, 0.5 * np.ones(n)]),
        np.column_stack([-0.5 * np.ones(n), -t]),
    ])

    normales = np.repeat(
        [[0, -1], [1, 0], [0, 1], [-1, 0]], n, axis=0
    )
    return points, normales


def calculer_divergence(points, normales, sigma, X, Y):

    g = np.zeros_like(X)

    for point, normale in zip(points, normales):
        dx = X - point[0]
        dy = Y - point[1]
        distance2 = dx**2 + dy**2

        gaussienne = np.exp(-distance2 / (2 * sigma**2))
        gaussienne /= 2 * np.pi * sigma**2

        produit = normale[0] * dx + normale[1] * dy
        g += produit * gaussienne

    return -g / sigma**2


x = np.linspace(-0.7, 0.7, 121)
y = np.linspace(-0.7, 0.7, 121)
X, Y = np.meshgrid(x, y)
h = x[1] - x[0]

points, normales = creer_carre()
generateur = np.random.default_rng(4)
direction_bruit = generateur.normal(size=points.shape)

sigma_sur_h = np.linspace(0.5, 6, 12)
niveaux_bruit = [0.5, 1, 2]  

sensibilites = {bruit: [] for bruit in niveaux_bruit}
etalements = []

for rapport in sigma_sur_h:
    sigma = rapport * h
    g_propre = calculer_divergence(points, normales, sigma, X, Y)

    seuil = 0.1 * np.max(np.abs(g_propre))
    etalements.append(np.mean(np.abs(g_propre) > seuil))

    for bruit in niveaux_bruit:
        points_bruites = points + bruit * h * direction_bruit
        g_bruite = calculer_divergence(
            points_bruites, normales, sigma, X, Y
        )
        erreur = np.mean(np.abs(g_bruite - g_propre))
        reference = np.mean(np.abs(g_propre))
        sensibilites[bruit].append(erreur / reference)


sensibilite_moyenne = np.mean([np.array(sensibilites[b]) for b in niveaux_bruit], axis=0)
ecart = np.abs(sensibilite_moyenne - np.array(etalements))

sigma_choisi = sigma_sur_h[np.argmin(ecart)]


fig, axes = plt.subplots(1, 2, figsize=(10, 4))

for bruit in niveaux_bruit:
    axes[0].plot(
        sigma_sur_h,
        100 * np.array(sensibilites[bruit]),
        marker="o",
        label=f"bruit = {bruit}h",
    )

axes[0].set_title("Sensibilite au bruit")
axes[0].set_xlabel(r"$\sigma/h$")
axes[0].set_ylabel("variation de la divergence (%)")
axes[0].grid(alpha=0.3)
axes[0].legend()
axes[0].axvline(sigma_choisi, color="black", linestyle="--")

axes[1].plot(
    sigma_sur_h,
    100 * np.array(etalements),
    marker="o",
    color="tab:orange",
)
axes[1].set_title("Etalement du champ")
axes[1].set_xlabel(r"$\sigma/h$")
axes[1].set_ylabel("proportion de la grille active (%)")
axes[1].grid(alpha=0.3)
axes[1].axvline(
    sigma_choisi,
    color="black",
    linestyle="--",
    label=rf"compromis : $\sigma={sigma_choisi:.1f}h$",
)
axes[1].legend()

fig.suptitle("Choix de sigma : stabilite contre etalement")
fig.tight_layout()
fig.savefig("choix_sigma_simple.png", dpi=200)

fig2, ax2 = plt.subplots(figsize=(6, 4))

ax2.plot(
    sigma_sur_h, 
    100 * ecart, 
    marker="o", 
    color="purple", 
    label="Écart (|Sensibilité - Étalement|)"
)

ax2.axvline(
    sigma_choisi, 
    color="black", 
    linestyle="--", 
    label=rf"Minimum optimal : $\sigma={sigma_choisi:.1f}h$"
)

ax2.set_title("Recherche automatique du compromis")
ax2.set_xlabel(r"$\sigma/h$")
ax2.set_ylabel("Écart entre les courbes (%)")
ax2.grid(alpha=0.3)
ax2.legend()
plt.tight_layout()
fig2.savefig("courbe_ecart_sigma.png", dpi=200)

plt.show()