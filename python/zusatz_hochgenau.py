"""Zusatz (nicht Teil der Aufgabe): Kontrollen mit 120 Stellen Genauigkeit.

Die Kennlinie ist die NIST-Formel E(T) für Typ K; T(E) wird durch Nullstellensuche bestimmt.
Fragen: 1) Ist die Stagnation von e_interp bei äquidistanten Stützstellen echt oder Rundung?
        2) Wie genau ist kappa aus cond() für riesige Konditionszahlen?
        3) Was sieht das Prüfraster nicht (Lambda am Rand)?
        4) Hält die gewählte Lösung (Tschebyscheff, n = 15) auch auf ganz [a, b]?
Dauert wenige Sekunden.
"""
import mpmath as mp

from interpolation import *

mp.mp.dps = 120
C = [mp.mpf(s) for s in (
    "-0.176004136860E-01 0.389212049750E-01 0.185587700320E-04 -0.994575928740E-07 "
    "0.318409457190E-09 -0.560728448890E-12 0.560750590590E-15 -0.320207200030E-18 "
    "0.971511471520E-22 -0.121047212750E-25").split()]
A0, A1, A2 = mp.mpf("0.118597600000E+00"), mp.mpf("-0.118343200000E-03"), mp.mpf("0.126968600000E+03")


def E_von_T(T):
    return mp.polyval(C[::-1], T) + A0 * mp.exp(A1 * (T - A2) ** 2)


def T_von_E(E, T_start):
    return mp.findroot(lambda T: E_von_T(T) - E, mp.mpf(float(T_start)))


E, T = lese_pruefraster()
E_teil, T_teil = E[::10], T[::10]       # 200 Punkte reichen hier
T_wahr = [T_von_E(mp.mpf(float(e)), t) for e, t in zip(E_teil, T_teil)]
print("Prüfraster gegen NIST-Formel: max. Abweichung "
      f"{max(abs(float(a - mp.mpf(float(b)))) for a, b in zip(T_wahr, T_teil)):.1e} °C\n")

# feines Gitter über ganz [a, b]
b = lese_tabelle("stuetzwerte_aequidistant.csv")[:, 2].max()
fein = np.linspace(0, b, 20001)

for datei in ("stuetzwerte_aequidistant.csv", "stuetzwerte_tschebyscheff.csv"):
    M = lese_tabelle(datei)
    print(datei)
    print("  n   e_interp 120 St.  e_interp double  kappa wahr  kappa double  Lambda [a,b]  Lambda Raster")
    for n in (2, 10, 15, 20, 25, 30, 35, 40):
        x, T_ex, _ = stuetzwerte(M, n)
        xm = [mp.mpf(float(v)) for v in x]
        y = [T_von_E(v, t) for v, t in zip(xm, T_ex)]
        w = [1 / mp.fprod(xm[j] - xm[k] for k in range(n + 1) if k != j) for j in range(n + 1)]

        e_hg = 0
        for e, t in zip(E_teil, T_wahr):
            z = mp.mpf(float(e))
            q = [w[j] / (z - xm[j]) for j in range(n + 1)]
            e_hg = max(e_hg, abs(mp.fsum(q[j] * y[j] for j in range(n + 1)) / mp.fsum(q) - t))

        wd = baryz_gewichte(x)
        e_d = np.max(np.abs(baryz_auswerten(x, T_ex, wd, E_teil) - T_teil))
        V = mp.matrix([[v ** k for k in range(n + 1)] for v in xm])
        kappa = mp.norm(V, mp.inf) * mp.norm(V ** -1, mp.inf)
        punkte = fein[np.min(np.abs(fein[:, None] - x), axis=1) > 1e-12]
        print(f"  {n:2d}  {float(e_hg):14.2e}  {e_d:15.2e}  {float(kappa):10.2e}  "
              f"{np.linalg.cond(vandermonde_matrix(x), np.inf):12.2e}  "
              f"{lebesgue_konstante(x, wd, punkte):12.3e}  {lebesgue_konstante(x, wd, E):13.3e}")
    print()

# Wachstum / Abfall pro Grad
Ra = experiment("stuetzwerte_aequidistant.csv", E, T)
Rt = experiment("stuetzwerte_tschebyscheff.csv", E, T)
g = NS >= 10
for name, R in (("äquidistant", Ra), ("Tschebyscheff", Rt)):
    print(f"{name}: e_interp fällt pro Grad um Faktor {10 ** -np.polyfit(NS[g], np.log10(R['e_interp'][g]), 1)[0]:.2f}, "
          f"kappa wächst um Faktor {10 ** np.polyfit(NS[g], np.log10(R['kappa'][g]), 1)[0]:.1f}")

# Lage des Interpolationsfehlers bei äquidistanten Stützstellen
x, T_ex, _ = stuetzwerte(lese_tabelle("stuetzwerte_aequidistant.csv"), 40)
f = np.abs(baryz_auswerten(x, T_ex, baryz_gewichte(x), E) - T)
print(f"äquidistant, n = 40, exakte Daten: Maximum bei E = {E[np.argmax(f)]:.3f} mV; "
      f"mittlere Hälfte: {f[(E > 5.16) & (E < 15.48)].max():.1e} °C")

# Lebesgue-Formeln der Vorlesung auf ganz [a, b]
ok = True
for n in NS:
    x, _, _ = stuetzwerte(lese_tabelle("stuetzwerte_tschebyscheff.csv"), n)
    L = lebesgue_konstante(x, baryz_gewichte(x), fein[np.min(np.abs(fein[:, None] - x), axis=1) > 1e-12])
    ok = ok and L <= 2 / np.pi * np.log(n + 1) + 1
print(f"Lambda_n(X_T) <= (2/pi) ln(n+1) + 1 für alle n auf [a, b]: {ok}")

# Gewählte Lösung auf ganz [a, b]: exakte Werte über Newton an der NIST-Formel
x, _, T_kal = stuetzwerte(lese_tabelle("stuetzwerte_tschebyscheff.csv"), 15)
wd = baryz_gewichte(x)
punkte = fein[np.min(np.abs(fein[:, None] - x), axis=1) > 1e-12]
T_fein = np.interp(punkte, E, T)     # Startwerte
T_fein = np.array([float(T_von_E(mp.mpf(float(e)), t)) for e, t in zip(punkte[::20], T_fein[::20])])
f = np.abs(baryz_auswerten(x, T_kal, wd, punkte[::20]) - T_fein)
print(f"Tschebyscheff, n = 15, baryzentrisch auf {len(T_fein)} Punkten in [0, {b:.3f}] mV: max. Fehler {f.max():.4f} °C")
