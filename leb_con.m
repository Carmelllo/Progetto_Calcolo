function L = leb_con(z, x)
% LEB_CON  Calcolo della costante di Lebesgue tramite formula baricentrica
%
% INPUT:
%   z : vettore riga dei nodi di interpolazione
%       z = [z_1, z_2, ..., z_n]
%
%   x : vettore colonna dei punti in cui si valuta la funzione di Lebesgue
%       (tipicamente una griglia fitta sull'intervallo di interesse)
%
% OUTPUT:
%   L : approssimazione della costante di Lebesgue,
%       ottenuta come massimo della funzione di Lebesgue sui punti x


% Numero di nodi di interpolazione
n = length(z);


% --------------------------------------------------
% Calcolo dei pesi baricentrici
% --------------------------------------------------
% Inizializza il vettore dei pesi come vettore riga di 1
w = ones(1, n);

% Ciclo su tutti i nodi
for i = 1:n
    % Calcolo del peso baricentrico w(i):
    %  w(i) = 1 / prod_{j != i} (z(i) - z(j))
    %
    % z([1:i-1, i+1:end]) seleziona tutti i nodi tranne z(i)
    % prod(...) calcola il prodotto di tutti gli elementi del vettore
    w(i) = 1 / prod(z(i) - z([1:i-1, i+1:end]));
end


% --------------------------------------------------
% Valutazione della funzione di Lebesgue
% --------------------------------------------------
% Inizializza il vettore che conterrà i valori della funzione di Lebesgue
% lambda(k) = Lambda(x(k))
lambda = zeros(length(x), 1);

% Ciclo su tutti i punti x
for k = 1:length(x)
    
    % Calcola le differenze tra il punto x(k) e tutti i nodi z
    % diff(i) = x(k) - z(i)
    diff = x(k) - z;
    
    % Controllo se x(k) coincide con uno dei nodi z
    % (in questo caso una differenza è zero)
    if any(diff == 0)
        
        % Se x(k) è un nodo, la funzione di Lebesgue vale esattamente 1
        lambda(k) = 1;
        
    else
        
        % Calcola i termini |w_i / (x(k) - z_i)| per tutti i nodi
        % "./" indica divisione elemento per elemento
        temp = abs(w ./ diff);
        
        % Formula baricentrica della funzione di Lebesgue:
        %
        % Lambda(x) = ( somma_i |w_i / (x - z_i)| )
        %             --------------------------------
        %             | somma_i (w_i / (x - z_i)) |
        lambda(k) = sum(temp) / abs(sum(w ./ diff));
        
    end
end


% --------------------------------------------------
% Costante di Lebesgue
% --------------------------------------------------
% La costante di Lebesgue è il massimo della funzione di Lebesgue
% sui punti x considerati
L = max(lambda);

end
