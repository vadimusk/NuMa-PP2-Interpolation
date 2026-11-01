# Programmierprojekt 2 – Wie heiß ist die Lötspitze?

Polynominterpolation einer Thermoelement-Kennlinie. Numerische Mathematik für Elektrotechnik, WS 2026/27.
Zwei gleichwertige Umsetzungen: MATLAB (R2026a, keine Toolboxen) und Python (NumPy, Matplotlib).

## Ausführen

**MATLAB** im Ordner `matlab/`:

```matlab
aufgabe1   % Vandermonde und baryzentrisch, äquidistant
aufgabe2   % dasselbe mit Tschebyscheff-Stützstellen
aufgabe3   % Fehlerzerlegung 3a, 3b, 3c und Zusatz zur Anmerkung
```

**Python** (einmalig `python3 -m venv .venv && .venv/bin/pip install -r requirements.txt`), im Ordner `python/`:

```bash
../.venv/bin/python aufgabe1.py
../.venv/bin/python aufgabe2.py
../.venv/bin/python aufgabe3.py
../.venv/bin/python zusatz_hochgenau.py   # optional, Kontrollen mit 120 Stellen
```

Die Plots landen in `plots/matlab/` bzw. `plots/python/`.

## Dateien

| MATLAB | Python | Inhalt |
|---|---|---|
| `aufgabe1.m` | `aufgabe1.py` | 1a, 1b, Abweichung für n = 20 |
| `aufgabe2.m` | `aufgabe2.py` | Tschebyscheff, gemeinsamer Plot |
| `aufgabe3.m` | `aufgabe3.py` | 3a, 3b, 3c, Zusatz zur Anmerkung |
| `experiment.m` | `interpolation.py` | Rechnet alle n = 2…40 für einen Stützstellensatz durch |
| `vandermonde_matrix.m`, `horner.m` | `interpolation.py` | Monombasis |
| `baryz_gewichte.m`, `baryz_auswerten.m` | `interpolation.py` | Baryzentrische Formel |
| `lebesgue_konstante.m` | `interpolation.py` | Λₙ(X) = max Σ\|ℓⱼ\| über dem Prüfraster |
| `lese_stuetzwerte.m`, `speichere_plot.m`, `farben.m`, `log_achse.m` | `interpolation.py` | Hilfsfunktionen |
| – | `zusatz_hochgenau.py` | Kontrollen mit NIST-Formel und 120 Stellen (nicht Teil der Aufgabe) |

## Ergebnisse

- **1a:** e_max fällt von 2,2 °C auf das Minimum 0,123 °C bei n = 14 und steigt ab n = 16 sprunghaft an. **Das Ziel 0,1 °C wird nicht erreicht.** κ∞(V) wächst exponentiell, bis etwa 10⁵⁸ bei n = 40.
- **1b:** Beide Verfahren sind bis n ≈ 18 praktisch gleich. Danach wächst der baryzentrische Fehler exponentiell auf 2,7·10⁶ °C bei n = 40. Für n = 20 liegt der größte Fehler am Rand: +10,98 °C bei E = 0,261 mV (in der mittleren Hälfte nur 0,006 °C).
- **2:** **Ziel erreicht ab n = 9** (0,047 °C). Ab n ≈ 15 Plateau bei 0,005–0,01 °C. Vandermonde weicht ab n = 21 um mehr als 1 % ab und liegt bei einzelnen großen n wieder über 0,1 °C (MATLAB: n = 32, 35, 40).
- **3a:** e_interp fällt bei Tschebyscheff um Faktor ≈ 1,7 pro Grad bis 5·10⁻⁹ °C, bei äquidistant nur um Faktor ≈ 1,08 bis 0,025 °C. Λ₄₀ ≈ 4,7·10⁹ gegenüber 2,9.
- **3b:** Bis n = 13 dominiert der Interpolationsfehler, ab n = 14 der Datenfehler (bei beiden Stützstellenwahlen). Bei Tschebyscheff bleibt e_daten zwischen 0,005 und 0,0095 °C: Der Fehler lässt sich nicht beliebig verkleinern.
- **3c:** Die Daten ändern sich relativ um 10⁻⁵, die Polynomwerte ebenfalls um ≈ 10⁻⁵, die Koeffizienten aber bis zu rund 100 %. e_koeff wächst um Faktor ≈ 3,8 pro Grad, κ∞ um ≈ 30, und überholt e_daten ab n = 24.
- **Fazit:** Tschebyscheff-Stützstellen, baryzentrische Auswertung, **n = 15** (e_max = 0,0047 °C). Das Minimum für das Ziel wäre n = 9; ab n ≈ 14 bringt ein höherer Grad nichts mehr, weil der Datenfehler dominiert.

## Hinweis zur Reproduzierbarkeit

Die baryzentrischen Werte e_max und e_daten stimmen zwischen MATLAB und Python auf mindestens 7 Stellen überein, e_interp auf mindestens 2 Stellen (am knappsten bei äquidistanten Stützstellen und großem n). Die Vandermonde-Kurven weichen ab n ≈ 21 voneinander ab, weil dort Rundungsdetails (Potenzen, LAPACK-Implementierung) das Ergebnis bestimmen. Das ist genau der Effekt aus der Anmerkung der Aufgabenstellung und kein Fehler im Code. Die genauen n-Werte in diesem Bereich sind deshalb nicht übertragbar.
