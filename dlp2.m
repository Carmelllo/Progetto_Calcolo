function dlp = dlp2(x, d)
% DLP2  Calcolo dei punti di Leja approssimati (algoritmo 2)
%       usando fattorizzazione LU della Vandermonde di Chebyshev
%
% INPUT:
%   x : vettore dei punti della mesh in [-1,1]
%   d : grado del polinomio interpolante
%
% OUTPUT:
%   dlp : vettore di d+1 nodi di Leja approssimati

    n = d + 1;
    M = length(x);

    V = zeros(M, n);
    for j = 1:n
        V(:, j) = cos((j-1) * acos(x));
    end

    [~, ~, p] = lu(V, 'vector');
    dlp = x(p(1:n));
end
