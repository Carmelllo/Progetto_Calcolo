function L = leb_con(z, x)
% LEB_CON  Calcolo della costante di Lebesgue
%
% INPUT:
%   z : vettore riga dei nodi di interpolazione (z_0,...,z_n)
%   x : vettore colonna dei punti in cui valutare la funzione di Lebesgue
%
% OUTPUT:
%   L : approssimazione della costante di Lebesgue

    n = length(z) - 1;
    m = length(x);

    lambda = zeros(m, 1);

    for i = 1:n+1
        li = ones(m, 1);

        for j = [1:i-1, i+1:n+1]
            li = li .* (x - z(j)) / (z(i) - z(j));
        end

        lambda = lambda + abs(li);
    end

    L = max(lambda);
end
