function L = lebesgue_konstante(x, w, E)
% LEBESGUE_KONSTANTE  Lambda_n(X) = max_E sum_j |l_j(E)| (Skript, Proposition 2.15),
% das Maximum wird über die Punkte E gebildet (hier das Prüfraster).
% Baryzentrisch: l_j(E) = (w_j / (E - x_j)) / sum_k (w_k / (E - x_k))

x = x(:);  w = w(:);  E = E(:);
Q = w.' ./ (E - x.');
L = max(sum(abs(Q), 2) ./ abs(sum(Q, 2)));
end
