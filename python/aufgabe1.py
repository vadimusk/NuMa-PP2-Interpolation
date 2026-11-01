"""Programmierprojekt 2 – Aufgabe 1: Vandermonde vs. baryzentrisch, äquidistante Stützstellen.

Daten: stuetzwerte_aequidistant.csv, y_i = T_kalibrator_C.
Fehlermaß: e_max = max_{E in Prüfraster} |p(E) - T_exakt(E)|
"""
from interpolation import *

E, T_E = lese_pruefraster()
R = experiment("stuetzwerte_aequidistant.csv", E, T_E)

# Aufgabe 1a: Vandermonde
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
ax1.plot(NS, R["e_vander"], "o-", color=VANDER, label="Vandermonde")
ziel_linie(ax1)
log_achse(ax1, "1a: $e_{max}(n)$, Vandermonde, äquidistant", r"$e_{max}$ in °C")
ax1.set_ylim(1e-2, 1e3)
ax1.legend(loc="upper left")
ax2.plot(NS, R["kappa"], "o-", color=NEUTRAL)
log_achse(ax2, "1a: Konditionszahl der Vandermonde-Matrix", r"$\kappa_\infty(V(X))$")
speichere_plot(fig, "a1a_vandermonde.png")

i_min = np.argmin(R["e_vander"])
print(f"1a) Vandermonde äquidistant: kleinster Fehler {R['e_vander'][i_min]:.3f} °C bei n = {NS[i_min]}")
print(f"    Ziel e_max < 0,1 °C erreicht: {'ja' if np.any(R['e_vander'] < 0.1) else 'nein'}")
print(f"    kappa: n = 2: {R['kappa'][0]:.1e},  n = 20: {R['kappa'][18]:.1e},  n = 40: {R['kappa'][38]:.1e}\n")

# Aufgabe 1b: baryzentrisch, beide Kurven in einem Plot
fig, ax = plt.subplots(figsize=(7.3, 4.4))
ax.plot(NS, R["e_vander"], "o-", color=VANDER, label="Vandermonde")
ax.plot(NS, R["e_baryz"], "s-", color=BARYZ, label="baryzentrisch")
ziel_linie(ax)
log_achse(ax, "1b: $e_{max}(n)$, äquidistante Stützstellen, Kalibratorwerte", r"$e_{max}$ in °C")
ax.set_ylim(1e-2, 1e7)
ax.legend(loc="upper left")
speichere_plot(fig, "a1b_vergleich.png")

i_min = np.argmin(R["e_baryz"])
print(f"1b) Baryzentrisch äquidistant: kleinster Fehler {R['e_baryz'][i_min]:.3f} °C bei n = {NS[i_min]}, "
      f"bei n = 40: {R['e_baryz'][-1]:.2e} °C")
print("    n   Vandermonde  baryzentrisch")
for n in (5, 10, 14, 15, 16, 18, 20, 25, 30, 40):
    print(f"    {n:2d}   {R['e_vander'][n - 2]:10.3e}   {R['e_baryz'][n - 2]:10.3e}")

# Aufgabe 1b: Abweichung für n = 20 über E
x, _, T_kal = stuetzwerte(lese_tabelle("stuetzwerte_aequidistant.csv"), 20)
abw = baryz_auswerten(x, T_kal, baryz_gewichte(x), E) - T_E

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
ax1.plot(E, abw, color=BARYZ, linewidth=1.2, label=r"$p(E) - T_{exakt}(E)$")
ax1.plot(x, 0 * x, "k.", label="Stützstellen")
ax1.set(xlabel="E in mV", ylabel="Abweichung in °C", title="1b: Abweichung für n = 20")
ax1.legend(loc="upper center")
ax2.semilogy(E, np.abs(abw), color=BARYZ, linewidth=1, label=r"$|p(E) - T_{exakt}(E)|$")
ziel_linie(ax2)
ax2.set(xlabel="E in mV", ylabel="Betrag in °C", title="1b: dasselbe logarithmisch")
ax2.legend(loc="upper center")
for ax in (ax1, ax2):
    ax.grid(True, color="#e1e0d9", linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)
speichere_plot(fig, "a1b_abweichung_n20.png")

i_max = np.argmax(np.abs(abw))
mitte = (E > 5.16) & (E < 15.48)    # mittlere Hälfte des Intervalls
print(f"\n    n = 20: größter Fehler {abw[i_max]:.2f} °C bei E = {E[i_max]:.3f} mV (T = {T_E[i_max]:.1f} °C), also am Rand")
print(f"    In der mittleren Hälfte [5,16; 15,48] mV ist der Fehler nur {np.abs(abw[mitte]).max():.4f} °C")
