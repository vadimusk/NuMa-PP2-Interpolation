function [x, T_exakt, T_kal] = lese_stuetzwerte(M, n)
% LESE_STUETZWERTE  Holt aus der eingelesenen Stützwert-Tabelle die n+1 Zeilen
% für den Polynomgrad n. Spalten von M: n, i, E_mV, T_exakt_C, T_kalibrator_C

zeilen = M(:, 1) == n;
x       = M(zeilen, 3);
T_exakt = M(zeilen, 4);
T_kal   = M(zeilen, 5);
end
