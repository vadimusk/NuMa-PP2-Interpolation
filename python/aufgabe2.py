"""Programmierprojekt 2 – Aufgabe 2: Tschebyscheff-Stützstellen (wieder Kalibratorwerte)."""
from interpolation import *

E, T_E = lese_pruefraster()
R = experiment("stuetzwerte_tschebyscheff.csv", E, T_E)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
ax1.plot(NS, R["e_vander"], "o-", color=VANDER, label="Vandermonde")
ax1.plot(NS, R["e_baryz"], "s-", color=BARYZ, label="baryzentrisch")
ziel_linie(ax1)
log_achse(ax1, "2: $e_{max}(n)$, Tschebyscheff, Kalibratorwerte", r"$e_{max}$ in °C")
ax1.legend(loc="upper right")
ax2.plot(NS, R["kappa"], "o-", color=NEUTRAL)
log_achse(ax2, "2: Konditionszahl der Vandermonde-Matrix", r"$\kappa_\infty(V(X))$")
speichere_plot(fig, "a2_tschebyscheff.png")

# Fragen
k0 = np.argmax(R["e_baryz"] < 0.1)
n0 = NS[k0]
print(f"Ziel e_max < 0,1 °C: baryzentrisch ab n = {n0} (e = {R['e_baryz'][k0]:.3f} °C)")
print(f"    baryzentrisch unter dem Ziel für alle n >= {n0}: {np.all(R['e_baryz'][k0:] < 0.1)}")

n_ueber = NS[(NS > n0) & (R["e_vander"] >= 0.1)]
print(f"Vandermonde ab n = {NS[np.argmax(R['e_vander'] < 0.1)]} unter dem Ziel, aber wieder darüber bei n = {n_ueber}")

# Ab wann unterscheiden sich die Verfahren? (relativer Unterschied der Fehler > 1 %)
rel = np.abs(R["e_vander"] - R["e_baryz"]) / R["e_baryz"]
print(f"Verfahren unterscheiden sich um mehr als 1 % ab n = {NS[np.argmax(rel > 0.01)]}")

i_min = np.argmin(R["e_baryz"])
print(f"Kleinster baryzentrischer Fehler: {R['e_baryz'][i_min]:.2e} °C bei n = {NS[i_min]}; "
      f"bei n = 40: {R['e_baryz'][-1]:.2e} °C")

print("\n    n   Vandermonde  baryzentrisch     kappa")
for n in (5, 8, 9, 10, 12, 15, 18, 20, 22, 24, 25, 30, 33, 36, 40):
    k = n - 2
    print(f"    {n:2d}   {R['e_vander'][k]:10.3e}   {R['e_baryz'][k]:10.3e}   {R['kappa'][k]:9.2e}")
