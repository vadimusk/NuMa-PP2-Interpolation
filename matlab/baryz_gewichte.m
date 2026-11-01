function w = baryz_gewichte(x)
% BARYZ_GEWICHTE  Baryzentrische Gewichte w_j = 1 / prod_{k ~= j} (x_j - x_k).
% Hängen nur von den Stützstellen ab, einmalig O(n^2).

x = x(:);
w = zeros(size(x));
for j = 1:numel(x)
    d = x(j) - x;
    d(j) = [];              % Faktor k = j auslassen
    w(j) = 1 / prod(d);
end
end
