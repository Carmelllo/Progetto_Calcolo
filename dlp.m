function dlp = dlp(x, d)
% DLP  Calcolo dei punti di Leja approssimati (algoritmo 1)
%
% INPUT:
%   x : vettore dei punti della mesh (equispaziati in [-1,1])
%   d : grado del polinomio interpolante
%
% OUTPUT:
%   dlp : vettore di d+1 nodi di Leja approssimati

    n = d + 1;
    dlp = zeros(n, 1);
    dlp(1) = x(1);

    for k = 2:n
        prod_val = ones(length(x), 1);
        for j = 1:k-1
            prod_val = prod_val .* abs(x - dlp(j));
        end

        [~, idx] = max(prod_val);
        dlp(k) = x(idx);
    end
end
