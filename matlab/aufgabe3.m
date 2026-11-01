%% Programmierprojekt 2 – Aufgabe 3: Fehlerzerlegung
% p_Verfahren - T_exakt = (p[T_exakt] - T_exakt)        Interpolationsfehler (3a)
%                       + (p[T_kal] - p[T_exakt])       Datenfehler          (3b)
%                       + (p_Verfahren - p[T_kal])      Koeffizientenfehler  (3c)
% Alle Maxima laufen über das Prüfraster. Im Ordner matlab/ ausführen.

clear; close all; clc;
c = farben();

P = readmatrix('../daten/pruefraster.csv', 'FileType', 'text', 'CommentStyle', '#', 'Delimiter', ' ');
E = P(:, 1);  T_E = P(:, 2);
Ra = experiment('stuetzwerte_aequidistant.csv', E, T_E);
Rt = experiment('stuetzwerte_tschebyscheff.csv', E, T_E);
n  = Ra.n;


%% Aufgabe 3a: Interpolationsfehler (exakte Daten, baryzentrisch)
L_aequi_formel = 2.^(n + 1) ./ (exp(1) * n .* log(n));    % Vorlesung, für große n
L_tsch_formel  = 2/pi * log(n + 1) + 1;                   % Vorlesung

figure('Name', 'Aufgabe 3a', 'Position', [100 100 1050 400]);
subplot(1, 2, 1);
semilogy(n, Ra.e_interp, 'o-', 'Color', c.aequi, 'LineWidth', 1.5, 'MarkerSize', 4); hold on;
semilogy(n, Rt.e_interp, 's-', 'Color', c.tscheb, 'LineWidth', 1.5, 'MarkerSize', 4);
legend('äquidistant', 'Tschebyscheff', 'Location', 'southwest');
log_achse('3a: Interpolationsfehler e_{interp}(n)', 'e_{interp} in °C');

subplot(1, 2, 2);
semilogy(n, Ra.lebesgue, 'o-', 'Color', c.aequi, 'LineWidth', 1.5, 'MarkerSize', 4); hold on;
semilogy(n, L_aequi_formel, '--', 'Color', c.neutral, 'LineWidth', 1);
semilogy(n, Rt.lebesgue, 's-', 'Color', c.tscheb, 'LineWidth', 1.5, 'MarkerSize', 4);
semilogy(n, L_tsch_formel, ':', 'Color', c.neutral, 'LineWidth', 1.5);
legend('\Lambda_n äquidistant', '2^{n+1} / (e n ln n)', '\Lambda_n Tschebyscheff', ...
    '(2/\pi) ln(n+1) + 1', 'Location', 'northwest');
log_achse('3a: Lebesgue-Konstante (auf dem Prüfraster)', '\Lambda_n(X)');
speichere_plot('a3a_interpolationsfehler.png');

fprintf('3a)  n   e_interp äqui  e_interp Tsch   Lambda äqui  Lambda Tsch\n');
for m = [2 5 10 15 20 25 30 35 40]
    k = m - 1;
    fprintf('    %2d    %10.2e     %10.2e    %9.2e    %9.2f\n', m, Ra.e_interp(k), Rt.e_interp(k), ...
        Ra.lebesgue(k), Rt.lebesgue(k));
end

% Lebesgue-Abschätzung: Der Bestapproximationsfehler ist unbekannt, aber mit Tschebyscheff gilt
%    e_interp / (1 + Lambda) <= min_q ||T - q|| <= e_interp
% Damit sieht man, wie viel die Schranke für äquidistante Stützstellen taugt.
beste_oben  = Rt.e_interp(end);
beste_unten = Rt.e_interp(end) / (1 + Rt.lebesgue(end));
fprintf('\n    n = 40: Bestapproximationsfehler zwischen %.1e und %.1e °C\n', beste_unten, beste_oben);
fprintf('    Schranke für äquidistant: zwischen %.1f und %.1f °C, tatsächlich %.3f °C\n', ...
    (1 + Ra.lebesgue(end)) * beste_unten, (1 + Ra.lebesgue(end)) * beste_oben, Ra.e_interp(end));


%% Aufgabe 3b: Datenfehler
figure('Name', 'Aufgabe 3b', 'Position', [100 100 1050 400]);
subplot(1, 2, 1);
plot_3b(n, Ra, c);
title('3b: äquidistante Stützstellen');
subplot(1, 2, 2);
plot_3b(n, Rt, c);
title('3b: Tschebyscheff-Stützstellen');
speichere_plot('a3b_datenfehler.png');

fprintf('\n3b) Datenfehler größer als Interpolationsfehler ab:\n');
fprintf('    äquidistant:   n = %d\n', n(find(Ra.e_daten > Ra.e_interp, 1)));
fprintf('    Tschebyscheff: n = %d\n', n(find(Rt.e_daten > Rt.e_interp, 1)));
fprintf('    Schranke eingehalten: äqui %d, Tsch %d\n', ...
    all(Ra.e_daten <= 0.005 * Ra.lebesgue), all(Rt.e_daten <= 0.005 * Rt.lebesgue));
fprintf('    e_daten Tschebyscheff für n >= 15: zwischen %.4f und %.4f °C\n', ...
    min(Rt.e_daten(14:end)), max(Rt.e_daten(14:end)));

% Kontrolle der Zerlegung (Dreiecksungleichung)
fprintf('    Dreiecksungl. e_baryz <= e_interp + e_daten: äqui %d, Tsch %d\n', ...
    all(Ra.e_baryz <= (Ra.e_interp + Ra.e_daten) * (1 + 1e-9)), ...
    all(Rt.e_baryz <= (Rt.e_interp + Rt.e_daten) * (1 + 1e-9)));


%% Aufgabe 3c: Koeffizientenfehler (nur Tschebyscheff)
daten_rel = 0.005 / 500;              % relative Änderung der Daten: 0,005 °C von bis zu 500 °C

figure('Name', 'Aufgabe 3c', 'Position', [100 100 1050 400]);
subplot(1, 2, 1);
semilogy(n, Rt.koeff_rel, 'o-', 'Color', c.vander, 'LineWidth', 1.5, 'MarkerSize', 4); hold on;
semilogy(n, Rt.e_daten / max(abs(T_E)), 's-', 'Color', c.daten, 'LineWidth', 1.5, 'MarkerSize', 4);
yline(daten_rel, '--', 'Color', c.neutral, 'LineWidth', 1);
legend('Koeffizienten: ||a~ - a|| / ||a||', 'Polynomwerte: e_{daten} / 500 °C', ...
    'Daten: ||y~ - y|| / ||y||', 'Location', 'southeast');
log_achse('3c: relative Änderung, Tschebyscheff', 'relative Änderung');
ylim([1e-8 1e2]);

subplot(1, 2, 2);
semilogy(n, Rt.kappa, 'o-', 'Color', c.neutral, 'LineWidth', 1.5, 'MarkerSize', 4); hold on;
semilogy(n, Rt.e_koeff, 's-', 'Color', c.vander, 'LineWidth', 1.5, 'MarkerSize', 4);
semilogy(n, Rt.e_daten, '^-', 'Color', c.daten, 'LineWidth', 1.5, 'MarkerSize', 4);
legend('\kappa_\infty(V(X))', 'e_{koeff} in °C', 'e_{daten} in °C', 'Location', 'northwest');
log_achse('3c: Koeffizientenfehler und Konditionszahl', 'Wert (logarithmisch)');
speichere_plot('a3c_koeffizientenfehler.png');

fprintf('\n3c)  n   ||a~-a||/||a||  e_daten     e_koeff      kappa\n');
for m = [5 10 15 20 22 24 25 30 33 36 40]
    k = m - 1;
    fprintf('    %2d    %10.2e  %9.2e  %9.2e  %9.2e\n', m, Rt.koeff_rel(k), Rt.e_daten(k), Rt.e_koeff(k), Rt.kappa(k));
end
fprintf('    Relative Datenänderung: %.1e\n', daten_rel);
fprintf('    e_koeff > e_daten ab n = %d\n', n(find(Rt.e_koeff > Rt.e_daten, 1)));

% Wachstum pro Grad n (n = 5..24, danach sättigt e_koeff)
g = n >= 5 & n <= 24;
s_kappa = polyfit(n(g), log10(Rt.kappa(g)), 1);
s_koeff = polyfit(n(g), log10(Rt.e_koeff(g)), 1);
fprintf('    Wachstum pro Grad: kappa Faktor %.1f, e_koeff Faktor %.1f\n', 10^s_kappa(1), 10^s_koeff(1));
fprintf('    e_koeff / (1e-16 * kappa) bei n = 20: %.1e\n', Rt.e_koeff(19) / (1e-16 * Rt.kappa(19)));

% p_Vandermonde interpoliert y + r (r = Residuum), also p_Vandermonde - p[y] = p[r]
%   ->  e_koeff <= Lambda_n * ||r||_inf
g = n >= 10;
fprintf('    Kontrolle e_koeff <= Lambda * Residuum (n >= 10): größtes Verhältnis %.2f\n', ...
    max(Rt.e_koeff(g) ./ (Rt.lebesgue(g) .* Rt.residuum(g))));


%% Zusatz zur Anmerkung: Vandermonde bei äquidistanten Stützstellen
% Dieselbe Rechnung einmal in mV und einmal mit x skaliert auf [-1, 1].
M = readmatrix('../daten/stuetzwerte_aequidistant.csv', 'FileType', 'text', 'CommentStyle', '#', 'Delimiter', ' ');
b = max(M(:, 3));
e_skaliert = zeros(size(n));
alt = warning('off', 'MATLAB:nearlySingularMatrix');
alt2 = warning('off', 'MATLAB:illConditionedMatrix');
for k = 1:numel(n)
    [x, ~, T_kal] = lese_stuetzwerte(M, n(k));
    a = vandermonde_matrix(2 * x / b - 1) \ T_kal;
    e_skaliert(k) = max(abs(horner(a, 2 * E / b - 1) - T_E));
end
warning(alt);  warning(alt2);

figure('Name', 'Anmerkung', 'Position', [100 100 700 420]);
semilogy(n, Ra.e_baryz, 's-', 'Color', c.baryz, 'LineWidth', 1.5, 'MarkerSize', 4); hold on;
semilogy(n, Ra.e_vander, 'o-', 'Color', c.vander, 'LineWidth', 1.5, 'MarkerSize', 4);
semilogy(n, e_skaliert, '^', 'Color', c.interp, 'LineWidth', 1, 'MarkerSize', 7);
semilogy(n, Ra.residuum, 'd-', 'Color', c.neutral, 'LineWidth', 1, 'MarkerSize', 3);
legend('baryzentrisch', 'Vandermonde in mV', 'Vandermonde skaliert auf [-1, 1]', ...
    'Residuum |p(x_i) - y_i| (Rechnung in mV)', 'Location', 'northwest');
log_achse('Anmerkung: äquidistante Stützstellen', 'Fehler in °C');
speichere_plot('anmerkung_aequidistant.png');

fprintf('\nAnmerkung)  n   baryz       Vander mV   Vander skal.  Residuum mV\n');
for m = [15 20 25 30 35 40]
    k = m - 1;
    fprintf('           %2d   %9.2e   %9.2e   %9.2e     %9.2e\n', m, Ra.e_baryz(k), Ra.e_vander(k), e_skaliert(k), Ra.residuum(k));
end


function plot_3b(n, R, c)
semilogy(n, R.e_interp, 'o-', 'Color', c.interp, 'LineWidth', 1.5, 'MarkerSize', 4); hold on;
semilogy(n, R.e_daten, 's-', 'Color', c.daten, 'LineWidth', 1.5, 'MarkerSize', 4);
semilogy(n, 0.005 * R.lebesgue, '--', 'Color', c.neutral, 'LineWidth', 1.2);
legend('e_{interp}', 'e_{daten}', '\Lambda_n \cdot 0,005 °C', 'Location', 'best');
log_achse('', 'Fehler in °C');
end
