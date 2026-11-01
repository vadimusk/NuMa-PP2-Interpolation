function p = horner(a, x)
% HORNER  Wertet p(x) = a_0 + a_1 x + ... + a_n x^n aus (Skript, Bemerkung 2.7).
% a(1) = a_0, ..., a(end) = a_n. x darf ein Vektor sein.
% Klammern von innen nach außen: a_0 + x(a_1 + x(a_2 + ... + x a_n)), Aufwand O(n).

p = a(end) * ones(size(x));
for k = numel(a)-1:-1:1
    p = a(k) + x .* p;
end
end
