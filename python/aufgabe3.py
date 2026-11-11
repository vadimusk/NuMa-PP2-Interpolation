"""Programmierprojekt 2 – Aufgabe 3: Fehlerzerlegung.

p_Verfahren - T_exakt = (p[T_exakt] - T_exakt)      Interpolationsfehler (3a)
                      + (p[T_kal] - p[T_exakt])     Datenfehler          (3b)
                      + (p_Verfahren - p[T_kal])    Koeffizientenfehler  (3c)
Alle Maxima laufen über das Prüfraster.
"""
from interpolation import *

E, T_E = lese_pruefraster()
Ra = experiment("stuetzwerte_aequidistant.csv", E, T_E)
Rt = experiment("stuetzwerte_tschebyscheff.csv", E, T_E)
n = NS

# Aufgabe 3a: Interpolationsfehler (exakte Daten, baryzentrisch)
L_aequi_formel = 2.0 ** (n + 1) / (np.e * n * np.log(n))    # Vorlesung, für große n
L_tsch_formel = 2 / np.pi * np.log(n + 1) + 1               # Vorlesung

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
ax1.plot(n, Ra["e_interp"], "o-", color=AEQUI, label="äquidistant")
ax1.plot(n, Rt["e_interp"], "s-", color=TSCHEB, label="Tschebyscheff")
log_achse(ax1, "3a: Interpolationsfehler $e_{interp}(n)$", r"$e_{interp}$ in °C")
ax1.legend(loc="lower left")
ax2.plot(n, Ra["lebesgue"], "o-", color=AEQUI, label=r"$\Lambda_n$ äquidistant")
ax2.plot(n, L_aequi_formel, "--", color=NEUTRAL, linewidth=1, label=r"$2^{n+1} / (e\, n \ln n)$")
ax2.plot(n, Rt["lebesgue"], "s-", color=TSCHEB, label=r"$\Lambda_n$ Tschebyscheff")
ax2.plot(n, L_tsch_formel, ":", color=NEUTRAL, linewidth=1.5, label=r"$(2/\pi) \ln(n+1) + 1$")
log_achse(ax2, "3a: Lebesgue-Konstante (auf dem Prüfraster)", r"$\Lambda_n(X)$")
ax2.legend(loc="upper left")
speichere_plot(fig, "a3a_interpolationsfehler.png")

print("3a)  n   e_interp äqui  e_interp Tsch   Lambda äqui  Lambda Tsch")
for m in (2, 5, 10, 15, 20, 25, 30, 35, 40):
    k = m - 2
    print(f"    {m:2d}    {Ra['e_interp'][k]:10.2e}     {Rt['e_interp'][k]:10.2e}    "
          f"{Ra['lebesgue'][k]:9.2e}    {Rt['lebesgue'][k]:9.2f}")

# Lebesgue-Abschätzung: Der Bestapproximationsfehler ist unbekannt, aber mit Tschebyscheff gilt
#    e_interp / (1 + Lambda) <= min_q ||T - q|| <= e_interp
beste_oben = Rt["e_interp"][-1]
beste_unten = Rt["e_interp"][-1] / (1 + Rt["lebesgue"][-1])
print(f"\n    n = 40: Bestapproximationsfehler zwischen {beste_unten:.1e} und {beste_oben:.1e} °C")
print(f"    Schranke für äquidistant: zwischen {(1 + Ra['lebesgue'][-1]) * beste_unten:.1f} und "
      f"{(1 + Ra['lebesgue'][-1]) * beste_oben:.1f} °C, tatsächlich {Ra['e_interp'][-1]:.3f} °C")

# Aufgabe 3b: Datenfehler
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
for ax, R, titel in zip(axes, (Ra, Rt), ("äquidistante ", "Tschebyscheff-")):
    ax.plot(n, R["e_interp"], "o-", color=INTERP, label=r"$e_{interp}$")
    ax.plot(n, R["e_daten"], "s-", color=DATEN_F, label=r"$e_{daten}$")
    ax.plot(n, 0.005 * R["lebesgue"], "--", color=NEUTRAL, label=r"$\Lambda_n \cdot 0{,}005$ °C")
    log_achse(ax, f"3b: {titel}Stützstellen")
    ax.legend(loc="best")
speichere_plot(fig, "a3b_datenfehler.png")

print("\n3b) Datenfehler größer als Interpolationsfehler ab:")
print(f"    äquidistant:   n = {n[np.argmax(Ra['e_daten'] > Ra['e_interp'])]}")
print(f"    Tschebyscheff: n = {n[np.argmax(Rt['e_daten'] > Rt['e_interp'])]}")
print(f"    Schranke eingehalten: äqui {np.all(Ra['e_daten'] <= 0.005 * Ra['lebesgue'])}, "
      f"Tsch {np.all(Rt['e_daten'] <= 0.005 * Rt['lebesgue'])}")
print(f"    e_daten Tschebyscheff für n >= 15: zwischen {Rt['e_daten'][13:].min():.4f} und {Rt['e_daten'][13:].max():.4f} °C")

# Kontrolle der Zerlegung (Dreiecksungleichung)
print("    Dreiecksungl. e_baryz <= e_interp + e_daten: "
      f"äqui {np.all(Ra['e_baryz'] <= (Ra['e_interp'] + Ra['e_daten']) * (1 + 1e-9))}, "
      f"Tsch {np.all(Rt['e_baryz'] <= (Rt['e_interp'] + Rt['e_daten']) * (1 + 1e-9))}")

# Aufgabe 3c: Koeffizientenfehler (nur Tschebyscheff)
daten_rel = 0.005 / 500     # relative Änderung der Daten: 0,005 °C von bis zu 500 °C

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
ax1.plot(n, Rt["koeff_rel"], "o-", color=VANDER, label=r"Koeffizienten: $\|\tilde a - a\| / \|a\|$")
ax1.plot(n, Rt["e_daten"] / np.abs(T_E).max(), "s-", color=DATEN_F, label=r"Polynomwerte: $e_{daten}$ / 500 °C")
ax1.axhline(daten_rel, linestyle="--", linewidth=1, color=NEUTRAL, label=r"Daten: $\|\tilde y - y\| / \|y\|$")
log_achse(ax1, "3c: relative Änderung, Tschebyscheff", "relative Änderung")
ax1.set_ylim(1e-8, 1e2)
ax1.legend(loc="lower right")
ax2.plot(n, Rt["kappa"], "o-", color=NEUTRAL, label=r"$\kappa_\infty(V(X))$")
ax2.plot(n, Rt["e_koeff"], "s-", color=VANDER, label=r"$e_{koeff}$ in °C")
ax2.plot(n, Rt["e_daten"], "^-", color=DATEN_F, label=r"$e_{daten}$ in °C")
log_achse(ax2, "3c: Koeffizientenfehler und Konditionszahl", "Wert (logarithmisch)")
ax2.legend(loc="upper left")
speichere_plot(fig, "a3c_koeffizientenfehler.png")

print("\n3c)  n   ||a~-a||/||a||  e_daten     e_koeff      kappa")
for m in (5, 10, 15, 20, 22, 24, 25, 30, 33, 36, 40):
    k = m - 2
    print(f"    {m:2d}    {Rt['koeff_rel'][k]:10.2e}  {Rt['e_daten'][k]:9.2e}  {Rt['e_koeff'][k]:9.2e}  {Rt['kappa'][k]:9.2e}")
print(f"    Relative Datenänderung: {daten_rel:.1e}")
print(f"    e_koeff > e_daten ab n = {n[np.argmax(Rt['e_koeff'] > Rt['e_daten'])]}")

# Wachstum pro Grad n (n = 5..24, danach sättigt e_koeff)
g = (n >= 5) & (n <= 24)
s_kappa = np.polyfit(n[g], np.log10(Rt["kappa"][g]), 1)[0]
s_koeff = np.polyfit(n[g], np.log10(Rt["e_koeff"][g]), 1)[0]
print(f"    Wachstum pro Grad: kappa Faktor {10 ** s_kappa:.1f}, e_koeff Faktor {10 ** s_koeff:.1f}")
print(f"    e_koeff / (1e-16 * kappa) bei n = 20: {Rt['e_koeff'][18] / (1e-16 * Rt['kappa'][18]):.1e}")

# p_Vandermonde interpoliert y + r (r = Residuum), also p_Vandermonde - p[y] = p[r]
#   ->  e_koeff <= Lambda_n * ||r||_inf
g = n >= 10
print("    Kontrolle e_koeff <= Lambda * Residuum (n >= 10): größtes Verhältnis "
      f"{np.max(Rt['e_koeff'][g] / (Rt['lebesgue'][g] * Rt['residuum'][g])):.2f}")
