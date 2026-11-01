"""Gemeinsame Funktionen für Programmierprojekt 2 (entspricht den .m-Dateien in matlab/)."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

WURZEL = Path(__file__).resolve().parent.parent
DATEN = WURZEL / "daten"
PLOTS = WURZEL / "plots" / "python"
NS = np.arange(2, 41)

# Feste Farben: jede Größe hat in allen Plots dieselbe Farbe (wie matlab/farben.m)
VANDER, BARYZ = "#2a78d6", "#eb6834"
INTERP, DATEN_F = "#1baf7a", "#6250d6"
AEQUI, TSCHEB = "#eda100", "#008300"
NEUTRAL = "#52514e"

plt.rcParams.update({"figure.facecolor": "white", "axes.edgecolor": "#c3c2b7",
                     "text.color": "#0b0b0b", "axes.labelcolor": "#52514e",
                     "xtick.color": "#52514e", "ytick.color": "#52514e",
                     "legend.frameon": False, "lines.linewidth": 1.5, "lines.markersize": 4})


def lese_pruefraster():
    M = np.loadtxt(DATEN / "pruefraster.csv", comments="#")
    return M[:, 0], M[:, 1]


def lese_tabelle(datei):
    """Spalten: n, i, E_mV, T_exakt_C, T_kalibrator_C"""
    return np.loadtxt(DATEN / datei, comments="#")


def stuetzwerte(M, n):
    z = M[:, 0] == n
    return M[z, 2], M[z, 3], M[z, 4]


def vandermonde_matrix(x):
    """V(X) wie im Skript (2.9): Zeile i = (1, x_i, ..., x_i^n)."""
    return x[:, None] ** np.arange(len(x))


def horner(a, x):
    """p(x) = a_0 + x(a_1 + x(a_2 + ... + x a_n)), a[0] = a_0 (Skript, Bemerkung 2.7)."""
    p = np.full_like(x, a[-1], dtype=float)
    for k in range(len(a) - 2, -1, -1):
        p = a[k] + x * p
    return p


def baryz_gewichte(x):
    """w_j = 1 / prod_{k != j} (x_j - x_k), einmalig O(n^2)."""
    w = np.zeros_like(x)
    for j in range(len(x)):
        w[j] = 1 / np.prod(np.delete(x[j] - x, j))
    return w


def baryz_auswerten(x, y, w, E):
    """p(E) = sum_j (w_j y_j / (E - x_j)) / sum_j (w_j / (E - x_j)), E != x_j."""
    Q = w / (E[:, None] - x)     # Zeile = Auswertungspunkt, Spalte = Stützstelle
    return (Q @ y) / Q.sum(axis=1)


def lebesgue_konstante(x, w, E):
    """Lambda_n(X) = max_E sum_j |l_j(E)| (Skript, Proposition 2.15), Maximum über E."""
    Q = w / (E[:, None] - x)
    return np.max(np.abs(Q).sum(axis=1) / np.abs(Q.sum(axis=1)))


def experiment(datei, E, T_E):
    """Rechnet für einen Stützstellensatz alle n = 2..40 durch (wie matlab/experiment.m)."""
    M = lese_tabelle(datei)
    R = {k: np.zeros(len(NS)) for k in ("e_vander", "e_baryz", "kappa", "e_interp", "e_daten",
                                        "lebesgue", "koeff_rel", "e_koeff", "residuum")}
    for k, n in enumerate(NS):
        x, T_ex, T_kal = stuetzwerte(M, n)

        V = vandermonde_matrix(x)
        a_kal = np.linalg.solve(V, T_kal)
        a_ex = np.linalg.solve(V, T_ex)
        p_vander = horner(a_kal, E)

        w = baryz_gewichte(x)
        p_kal = baryz_auswerten(x, T_kal, w, E)
        p_ex = baryz_auswerten(x, T_ex, w, E)

        R["e_vander"][k] = np.max(np.abs(p_vander - T_E))
        R["e_baryz"][k] = np.max(np.abs(p_kal - T_E))
        R["kappa"][k] = np.linalg.cond(V, np.inf)
        R["e_interp"][k] = np.max(np.abs(p_ex - T_E))
        R["e_daten"][k] = np.max(np.abs(p_kal - p_ex))
        R["lebesgue"][k] = lebesgue_konstante(x, w, E)
        R["koeff_rel"][k] = np.linalg.norm(a_kal - a_ex, np.inf) / np.linalg.norm(a_ex, np.inf)
        R["e_koeff"][k] = np.max(np.abs(p_vander - p_kal))
        R["residuum"][k] = np.max(np.abs(horner(a_kal, x) - T_kal))
    return R


def speichere_plot(fig, name):
    PLOTS.mkdir(parents=True, exist_ok=True)
    fig.savefig(PLOTS / name, dpi=150, bbox_inches="tight")
    plt.close(fig)


def log_achse(ax, titel, ylabel="Fehler in °C"):
    ax.set_yscale("log")
    ax.set_xlabel("Polynomgrad n")
    ax.set_ylabel(ylabel)
    ax.set_title(titel)
    ax.set_xlim(2, 40)
    ax.grid(True, which="major", color="#e1e0d9", linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)


def ziel_linie(ax):
    ax.axhline(0.1, linestyle="--", linewidth=1, color=NEUTRAL, label="Ziel 0,1 °C")
