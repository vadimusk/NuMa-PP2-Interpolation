%% Programmierprojekt 2 – Aufgabe 2: Tschebyscheff-Stützstellen
% Wie Aufgabe 1, aber mit stuetzwerte_tschebyscheff.csv (wieder Kalibratorwerte).
% Im Ordner matlab/ ausführen.

clear; close all; clc;
c = farben();

P = readmatrix('../daten/pruefraster.csv', 'FileType', 'text', 'CommentStyle', '#', 'Delimiter', ' ');
E = P(:, 1);  T_E = P(:, 2);
R = experiment('stuetzwerte_tschebyscheff.csv', E, T_E);


%% Gemeinsamer Plot
figure('Name', 'Aufgabe 2', 'Position', [100 100 1050 400]);
subplot(1, 2, 1);
semilogy(R.n, R.e_vander, 'o-', 'Color', c.vander, 'LineWidth', 1.5, 'MarkerSize', 4); hold on;
semilogy(R.n, R.e_baryz, 's-', 'Color', c.baryz, 'LineWidth', 1.5, 'MarkerSize', 4);
yline(0.1, '--', 'Color', c.neutral, 'LineWidth', 1);
legend('Vandermonde', 'baryzentrisch', 'Ziel 0,1 °C', 'Location', 'northeast');
log_achse('2: e_{max}(n), Tschebyscheff-Stützstellen, Kalibratorwerte', 'e_{max} in °C');

subplot(1, 2, 2);
semilogy(R.n, R.kappa, 'o-', 'Color', c.neutral, 'LineWidth', 1.5, 'MarkerSize', 4);
log_achse('2: Konditionszahl der Vandermonde-Matrix', '\kappa_\infty(V(X))');
speichere_plot('a2_tschebyscheff.png');


%% Fragen
n0 = R.n(find(R.e_baryz < 0.1, 1));
fprintf('Ziel e_max < 0,1 °C: baryzentrisch ab n = %d (e = %.3f °C)\n', n0, R.e_baryz(n0 - 1));
fprintf('    baryzentrisch unter dem Ziel für alle n >= %d: %d\n', n0, all(R.e_baryz(n0-1:end) < 0.1));

n_ueber = R.n(R.n > n0 & R.e_vander >= 0.1);
fprintf('Vandermonde ab n = %d unter dem Ziel, aber wieder darüber bei n =%s\n', ...
    R.n(find(R.e_vander < 0.1, 1)), sprintf(' %d', n_ueber));

% Ab wann unterscheiden sich die Verfahren? (relativer Unterschied der Fehler > 1 %)
rel = abs(R.e_vander - R.e_baryz) ./ R.e_baryz;
fprintf('Verfahren unterscheiden sich um mehr als 1 %% ab n = %d\n', R.n(find(rel > 0.01, 1)));

[~, i_min] = min(R.e_baryz);
fprintf('Kleinster baryzentrischer Fehler: %.2e °C bei n = %d; bei n = 40: %.2e °C\n', ...
    R.e_baryz(i_min), R.n(i_min), R.e_baryz(end));

fprintf('\n    n   Vandermonde  baryzentrisch     kappa\n');
for n = [5 8 9 10 12 15 18 20 22 24 25 30 33 36 40]
    k = n - 1;
    fprintf('    %2d   %10.3e   %10.3e   %9.2e\n', n, R.e_vander(k), R.e_baryz(k), R.kappa(k));
end
