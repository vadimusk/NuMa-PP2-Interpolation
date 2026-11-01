function R = experiment(datei, E, T_E)
% EXPERIMENT  Rechnet für einen Stützstellensatz alle n = 2..40 durch.
% Alle Maxima laufen über das Prüfraster E mit den exakten Werten T_E.
% Die Ergebnisse stehen als Spaltenvektoren (Länge 39) in der Struktur R.

M  = readmatrix(fullfile('..', 'daten', datei), 'FileType', 'text', ...
                'CommentStyle', '#', 'Delimiter', ' ');
ns = (2:40)';

R.n = ns;
R.e_vander = zeros(size(ns));   % |p_Vandermonde[T_kal] - T_exakt|
R.e_baryz  = zeros(size(ns));   % |p_baryz[T_kal] - T_exakt|
R.kappa    = zeros(size(ns));   % kappa_inf(V(X))
R.e_interp = zeros(size(ns));   % |p[T_exakt] - T_exakt|         (3a)
R.e_daten  = zeros(size(ns));   % |p[T_kal] - p[T_exakt]|        (3b)
R.lebesgue = zeros(size(ns));   % Lambda_n(X)
R.koeff_rel = zeros(size(ns));  % ||a_kal - a|| / ||a||          (3c)
R.e_koeff  = zeros(size(ns));   % |p_Vandermonde - p_baryz|      (3c)
R.residuum = zeros(size(ns));   % max_i |p_Vandermonde(x_i) - y_i|  (Anmerkung)

% Bei großem n warnt MATLAB zu Recht vor der fast singulären Matrix - das wollen wir
% hier ja gerade untersuchen, also Warnung für diese Funktion ausschalten.
alt = warning('off', 'MATLAB:nearlySingularMatrix');
alt2 = warning('off', 'MATLAB:illConditionedMatrix');

for k = 1:numel(ns)
    [x, T_ex, T_kal] = lese_stuetzwerte(M, ns(k));

    % Monombasis
    V = vandermonde_matrix(x);
    a_kal = V \ T_kal;
    a_ex  = V \ T_ex;
    p_vander = horner(a_kal, E);

    % Baryzentrisch
    w = baryz_gewichte(x);
    p_kal = baryz_auswerten(x, T_kal, w, E);
    p_ex  = baryz_auswerten(x, T_ex,  w, E);

    R.e_vander(k) = max(abs(p_vander - T_E));
    R.e_baryz(k)  = max(abs(p_kal - T_E));
    R.kappa(k)    = cond(V, inf);
    R.e_interp(k) = max(abs(p_ex - T_E));
    R.e_daten(k)  = max(abs(p_kal - p_ex));
    R.lebesgue(k) = lebesgue_konstante(x, w, E);
    R.koeff_rel(k) = norm(a_kal - a_ex, inf) / norm(a_ex, inf);
    R.e_koeff(k)  = max(abs(p_vander - p_kal));
    R.residuum(k) = max(abs(horner(a_kal, x) - T_kal));
end

warning(alt);  warning(alt2);
end
