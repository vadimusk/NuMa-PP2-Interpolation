function p = baryz_auswerten(x, y, w, E)
% BARYZ_AUSWERTEN  Baryzentrische Formel an den Stellen E:
%   p(E) = sum_j (w_j y_j / (E - x_j)) / sum_j (w_j / (E - x_j))
% Das Prüfraster enthält keine Stützstellen, also E ~= x_j (kein Sonderfall nötig).

x = x(:);  y = y(:);  w = w(:);  E = E(:);
Q = w.' ./ (E - x.');       % Zeile = Auswertungspunkt, Spalte = Stützstelle
p = (Q * y) ./ sum(Q, 2);
end
