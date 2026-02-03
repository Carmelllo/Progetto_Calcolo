function L = leb_con(z, x)
% LEB_CON  Calcolo della costante di Lebesgue tramite formula baricentrica
%
% INPUT:
%   z : vettore riga dei nodi di interpolazione%
%   x : vettore colonna dei punti in cui si valuta la funzione di Lebesgue
%
% OUTPUT:
%   L : approssimazione della costante di Lebesgue,
%       ottenuta come massimo della funzione di Lebesgue sui punti x


n = length(z);
w = ones(1, n);

for i = 1:n
    w(i) = 1 / prod(z(i) - z([1:i-1, i+1:end]));
end
lambda = zeros(length(x), 1);

for k = 1:length(x)
    diff = x(k) - z;    %vettore
    if any(diff == 0)
        lambda(k) = 1; 
    else
        temp = abs(w ./ diff);  
        lambda(k) = sum(temp) / abs(sum(w ./ diff)); 
    end
end

L = max(lambda);

end
