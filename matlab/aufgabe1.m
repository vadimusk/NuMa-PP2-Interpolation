%% Programmierprojekt 2 – Aufgabe 1: Vandermonde vs. baryzentrisch, äquidistante Stützstellen
% Daten: Stützstellen aus stuetzwerte_aequidistant.csv, y_i = T_kalibrator_C.
% Fehlermaß: e_max = max_{E in Prüfraster} |p(E) - T_exakt(E)|
% Im Ordner matlab/ ausführen.

clear; close all; clc;
c = farben();

P = readmatrix('../daten/pruefraster.csv', 'FileType', 'text', 'CommentStyle', '#', 'Delimiter', ' ');
E = P(:, 1);  T_E = P(:, 2);
R = experiment('stuetzwerte_aequidistant.csv', E, T_E);


%% Aufgabe 1a: Vandermonde
figure('Name', 'Aufgabe 1a', 'Position', [100 100 1050 400]);
subplot(1, 2, 1);
semilogy(R.n, R.e_vander, 'o-', 'Color', c.vander, 'LineWidth', 1.5, 'MarkerSize', 4); hold on;
yline(0.1, '--', 'Color', c.neutral, 'LineWidth', 1);
legend('Vandermonde', 'Ziel 0,1 °C', 'Location', 'northwest');
log_achse('1a: e_{max}(n), Vandermonde, äquidistant', 'e_{max} in °C');
ylim([1e-2 1e3]);

subplot(1, 2, 2);
semilogy(R.n, R.kappa, 'o-', 'Color', c.neutral, 'LineWidth', 1.5, 'MarkerSize', 4);
log_achse('1a: Konditionszahl der Vandermonde-Matrix', '\kappa_\infty(V(X))');
speichere_plot('a1a_vandermonde.png');

[e_min, i_min] = min(R.e_vander);
fprintf('1a) Vandermonde äquidistant: kleinster Fehler %.3f °C bei n = %d\n', e_min, R.n(i_min));
fprintf('    Ziel e_max < 0,1 °C erreicht: %s\n', ja_nein(any(R.e_vander < 0.1)));
fprintf('    kappa: n = 2: %.1e,  n = 20: %.1e,  n = 40: %.1e\n\n', R.kappa(1), R.kappa(19), R.kappa(39));


%% Aufgabe 1b: baryzentrisch, beide Kurven in einem Plot
figure('Name', 'Aufgabe 1b', 'Position', [100 100 700 420]);
semilogy(R.n, R.e_vander, 'o-', 'Color', c.vander, 'LineWidth', 1.5, 'MarkerSize', 4); hold on;
semilogy(R.n, R.e_baryz, 's-', 'Color', c.baryz, 'LineWidth', 1.5, 'MarkerSize', 4);
yline(0.1, '--', 'Color', c.neutral, 'LineWidth', 1);
legend('Vandermonde', 'baryzentrisch', 'Ziel 0,1 °C', 'Location', 'northwest');
log_achse('1b: e_{max}(n), äquidistante Stützstellen, Kalibratorwerte', 'e_{max} in °C');
ylim([1e-2 1e7]);
speichere_plot('a1b_vergleich.png');

[e_min, i_min] = min(R.e_baryz);
fprintf('1b) Baryzentrisch äquidistant: kleinster Fehler %.3f °C bei n = %d, bei n = 40: %.2e °C\n', ...
    e_min, R.n(i_min), R.e_baryz(end));
fprintf('    n   Vandermonde  baryzentrisch\n');
for n = [5 10 14 15 16 18 20 25 30 40]
    k = n - 1;
    fprintf('    %2d   %10.3e   %10.3e\n', n, R.e_vander(k), R.e_baryz(k));
end


%% Aufgabe 1b: Abweichung für n = 20 über E
M = readmatrix('../daten/stuetzwerte_aequidistant.csv', 'FileType', 'text', 'CommentStyle', '#', 'Delimiter', ' ');
[x, ~, T_kal] = lese_stuetzwerte(M, 20);
w = baryz_gewichte(x);
abw = baryz_auswerten(x, T_kal, w, E) - T_E;

figure('Name', 'Aufgabe 1b - n = 20', 'Position', [100 100 1050 400]);
subplot(1, 2, 1);
plot(E, abw, 'Color', c.baryz, 'LineWidth', 1.2); hold on;
plot(x, zeros(size(x)), 'k.', 'MarkerSize', 9);
legend('p(E) - T_{exakt}(E)', 'Stützstellen', 'Location', 'north');
grid on;  box off;
xlim([0 20.644]);
xlabel('E in mV');  ylabel('Abweichung in °C');
title('1b: Abweichung für n = 20');

subplot(1, 2, 2);
semilogy(E, abs(abw), 'Color', c.baryz, 'LineWidth', 1); hold on;
yline(0.1, '--', 'Color', c.neutral, 'LineWidth', 1);
legend('|p(E) - T_{exakt}(E)|', 'Ziel 0,1 °C', 'Location', 'north');
grid on;  box off;
set(gca, 'YMinorGrid', 'off');
xlim([0 20.644]);
xlabel('E in mV');  ylabel('Betrag in °C');
title('1b: dasselbe logarithmisch');
speichere_plot('a1b_abweichung_n20.png');

[e_max, i_max] = max(abs(abw));
mitte = E > 5.16 & E < 15.48;       % mittlere Hälfte des Intervalls
fprintf('\n    n = 20: größter Fehler %.2f °C bei E = %.3f mV (T = %.1f °C), also am Rand\n', ...
    abw(i_max), E(i_max), T_E(i_max));
fprintf('    In der mittleren Hälfte [5,16; 15,48] mV ist der Fehler nur %.4f °C\n', max(abs(abw(mitte))));


function s = ja_nein(b)
if b, s = 'ja'; else, s = 'nein'; end
end
