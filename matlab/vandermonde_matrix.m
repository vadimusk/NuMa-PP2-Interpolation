function V = vandermonde_matrix(x)
% VANDERMONDE_MATRIX  V(X) wie im Skript (2.9): Zeile i = (1, x_i, x_i^2, ..., x_i^n).
% Achtung: MATLABs eigenes vander() hat die Spalten in umgekehrter Reihenfolge.

x = x(:);
n = numel(x) - 1;
V = x .^ (0:n);
end
