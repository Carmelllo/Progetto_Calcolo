# Guida Completa al Progetto: Spiegazione Riga per Riga e Teoria

Questo documento è progettato per prepararti all'esame orale. Contiene:
1.  **Analisi delle Richieste:** Come il codice soddisfa i punti del PDF "Prog1".
2.  **Spiegazione Riga per Riga:** Analisi dettagliata di ogni singolo comando MATLAB dei 4 file `.m`.
3.  **Teoria e Grafici:** Interpretazione profonda dei risultati.

---

## PARTE 1: Risoluzione dei Problemi del PDF (Prog1)

Il PDF richiede di sviluppare un progetto per l'interpolazione polinomiale. Ecco come ogni richiesta è stata implementata:

### 1. "Utilizzare punti di Leja approssimati estratti da una mesh fitta"
**Richiesta:** Non calcolare i veri punti di Leja (complessi), ma sceglierli da un sottoinsieme discreto $K \subset [-1, 1]$.
**Implementazione:**
*   File: `main_experiment.m`, `dlp.m`, `dlp2.m`.
*   Codice: Creiamo `x_leja` con 100.000 punti equispaziati (la "mesh fitta"). Gli algoritmi `dlp` e `dlp2` estraggono i punti da questo insieme.

### 2. "Interpolante su base di Chebyshev di primo tipo"
**Richiesta:** Usare la formula $V_{ij} = \cos((j-1)\arccos(z_i))$.
**Implementazione:**
*   File: `dlp2.m` e `main_experiment.m`.
*   Codice: Non usiamo basi standard ($1, x, x^2$). Costruiamo esplicitamente la matrice $V$ usando `cos` e `acos` riga per riga nei cicli `for`.

### 3. "Risolvere il sistema lineare... senza usare polyfit"
**Richiesta:** Trovare i coefficienti $c$ risolvendo $Vc = f(z)$.
**Implementazione:**
*   File: `main_experiment.m`.
*   Codice: Usiamo l'operatore **backslash** (`\`) di MATLAB: `coeff = V \ y`. Questo è il metodo standard e stabile per sistemi lineari, evitando l'inversa e `polyfit`.

### 4. "Calcolare la costante di Lebesgue"
**Richiesta:** Valutare la stabilità $L = \max \sum |l_i(x)|$.
**Implementazione:**
*   File: `leb_con.m`.
*   Codice: Implementa direttamente la formula della sommatoria dei valori assoluti dei polinomi di Lagrange e ne calcola il massimo.

---

## PARTE 2: Analisi del Codice Riga per Riga

### File 1: `dlp.m` (Algoritmo Iterativo / Greedy)

Questo file implementa la definizione classica: "Scegli il prossimo punto che massimizza il prodotto delle distanze dai precedenti".

```matlab
function dlp = dlp(x, d)
```
**Spiegazione:** Definisce la funzione. Accetta la mesh `x` (i candidati) e il grado `d`.

```matlab
    n = d + 1;
```
**Spiegazione:** Se interpoliamo un polinomio di grado $d$, abbiamo bisogno di $d+1$ nodi. Esempio: grado 1 (retta) $\to$ 2 punti.

```matlab
    dlp = zeros(n, 1);
```
**Spiegazione:** **Preallocazione**. Crea un vettore colonna vuoto. Fondamentale per la velocità in MATLAB (evita di ingrandire il vettore a ogni ciclo).

```matlab
    dlp(1) = x(1);
```
**Spiegazione:** **Inizializzazione**. Il primo punto di Leja è arbitrario, ma solitamente si prende quello col modulo massimo. Qui prendiamo il primo elemento della mesh (che è -1).

```matlab
    for k = 2:n
```
**Spiegazione:** Ciclo principale. Inizia a cercare dal 2° punto fino all'ultimo ($n$).

```matlab
        prod_val = ones(length(x), 1);
```
**Spiegazione:** Crea un vettore di "1" lungo quanto la mesh (100.000 elementi). Questo vettore accumulerà il prodotto delle distanze. Si parte da 1 perché è l'elemento neutro della moltiplicazione.

```matlab
        for j = 1:k-1
```
**Spiegazione:** Ciclo interno. Scorre tutti i nodi di Leja *già trovati* finora.

```matlab
            prod_val = prod_val .* abs(x - dlp(j));
```
**Spiegazione:** **Riga Chiave**.
*   `x - dlp(j)`: Calcola la distanza tra *ogni* punto della mesh e il nodo $j$-esimo.
*   `abs(...)`: Valore assoluto (distanza positiva).
*   `.*`: **Moltiplicazione elemento-per-elemento**. Aggiorna il prodotto cumulativo per ogni candidato.
*   *Matematica:* Calcola $\prod_{j=0}^{k-1} |x - \xi_j|$.

```matlab
        end
        [~, idx] = max(prod_val);
```
**Spiegazione:**
*   `max(...)`: Trova il valore massimo nel vettore dei prodotti.
*   `[~, idx]`: Ci interessa solo **dove** si trova il massimo (l'indice `idx`), non quanto vale.

```matlab
        dlp(k) = x(idx);
```
**Spiegazione:** Aggiunge il punto trovato alla lista dei nodi scelti.

```matlab
    end
end
```

---

### File 2: `dlp2.m` (Algoritmo Matriciale / Algebra Lineare)

Questo file ottiene lo stesso risultato usando la fattorizzazione LU sulla matrice di Vandermonde. È più veloce per gradi alti.

```matlab
function dlp = dlp2(x, d)
    n = d + 1;
    M = length(x);
```
**Spiegazione:** Salva in `M` la grandezza della mesh (es. 100.000).

```matlab
    V = zeros(M, n);
```
**Spiegazione:** Prealloca la matrice $M \times n$. Molto grande in altezza, stretta in larghezza.

```matlab
    for j = 1:n
        V(:, j) = cos((j-1) * acos(x));
    end
```
**Spiegazione:** **Costruzione Matrice**.
*   `V(:, j)`: Accede all'intera colonna $j$-esima.
*   `cos((j-1)*acos(x))`: Applica i polinomi di Chebyshev. Se $j=1$, grado 0. Se $j=2$, grado 1, ecc.
*   Questa operazione è "vettorializzata": calcola il coseno per tutti i 100.000 punti in un colpo solo.

```matlab
    [~, ~, p] = lu(V, 'vector');
```
**Spiegazione:** **Fattorizzazione LU con Pivoting**.
*   Il comando `lu` scompone la matrice.
*   Il parametro `'vector'` chiede di restituire `p` come vettore di permutazione.
*   *Perché funziona?* L'algoritmo LU scambia le righe per mettere gli elementi "migliori" (più grandi) sulla diagonale per stabilità numerica. Le righe scelte per prime corrispondono ai punti di Leja.

```matlab
    dlp = x(p(1:n));
```
**Spiegazione:** Prende i primi $n$ indici dal vettore `p` e estrae i punti corrispondenti dalla mesh.

---

### File 3: `leb_con.m` (Costante di Lebesgue)

Calcola la stabilità dei nodi scelti.

```matlab
function L = leb_con(z, x)
    n = length(z) - 1;
    m = length(x);
    lambda = zeros(m, 1);
```
**Spiegazione:** Prepara il vettore `lambda` che accumulerà la somma.

```matlab
    for i = 1:n+1
        li = ones(m, 1);
```
**Spiegazione:** Inizia il calcolo dell'$i$-esimo polinomio di base di Lagrange $l_i(x)$.

```matlab
        for j = [1:i-1, i+1:n+1]
```
**Spiegazione:** Ciclo che scorre tutti gli indici $j$ **tranne** $i$. Serve per saltare il termine $(x-x_i)$ nel prodotto.

```matlab
            li = li .* (x - z(j)) / (z(i) - z(j));
```
**Spiegazione:** **Formula di Lagrange**.
*   Moltiplica accumulando il prodotto.
*   Nota: `z(i)` è uno scalare, `x` è un vettore. MATLAB gestisce automaticamente questa operazione ("broadcasting").

```matlab
        end
        lambda = lambda + abs(li);
```
**Spiegazione:** Somma il valore assoluto del polinomio appena calcolato alla funzione totale $\lambda(x)$.

```matlab
    end
    L = max(lambda);
```
**Spiegazione:** **Costante di Lebesgue**. È il valore massimo raggiunto dalla funzione somma sull'intervallo.

---

### File 4: `main_experiment.m` (Script Principale)

```matlab
clearvars; close all; clc
```
**Spiegazione:** Resetta l'ambiente di lavoro.

```matlab
M_leja = 100000;
x_leja = linspace(-1, 1, M_leja)';
```
**Spiegazione:** Crea la **Mesh Fitta** richiesta dal punto 1 del PDF. L'apice `'` traspone per avere un vettore colonna.

```matlab
M_eval = 2001;
x_eval = linspace(-1, 1, M_eval)';
```
**Spiegazione:** Crea una griglia di punti "test" per misurare l'errore e disegnare i grafici.

```matlab
degrees = 1:50;
f_test = @(x) 1./(x - 1.3);
```
**Spiegazione:**
*   Testiamo gradi da 1 a 50.
*   **Funzione anonima:** Definiamo $f(x) = \frac{1}{x-1.3}$. Abbiamo scelto questa funzione perché ha una singolarità vicino a 1, perfetta per testare il fenomeno di Runge.

```matlab
% Ciclo Principale
for ii = 1:length(degrees)
    d = degrees(ii);
    n = d + 1;
    
    % --- Algoritmo 1 ---
    tic;                 % Avvia cronometro
    z1 = dlp(x_leja, d); % Chiama funzione
    time_dlp(ii) = toc;  % Ferma cronometro e salva tempo
    L_dlp(ii) = leb_con(z1', x_eval); % Calcola Lebesgue
    
    % --- Algoritmo 2 ---
    ... (codice analogo per dlp2) ...
    
    % --- Interpolazione (Punti B e C del PDF) ---
    V_leja = zeros(n, n);
    V_eval = zeros(M_eval, n);
    for j = 1:n
        V_leja(:, j) = cos((j-1) * acos(z2));     % Matrice quadrata per coeff
        V_eval(:, j) = cos((j-1) * acos(x_eval)); % Matrice rettangolare per valutazione
    end
    
    coeff_leja = V_leja \ f_test(z2); % Risolve sistema lineare (no inversa!)
    p_leja = V_eval * coeff_leja;     % Calcola valori polinomio
    
    err_leja(ii) = max(abs(f_test(x_eval) - p_leja)); % Errore massimo
    
    % --- Confronto Equispaziati (Punto D del PDF) ---
    zeq = linspace(-1, 1, n)'; % Nodi equispaziati
    ... (stessa procedura di sopra) ...
end
```

### Codice Grafici (Visualizzazione)

```matlab
figure('Name', '...');
semilogy(degrees, L_dlp, '-b');
```
**Spiegazione:** `semilogy` crea un grafico con asse Y logaritmico. Essenziale per visualizzare l'errore che varia da $10^{-1}$ a $10^{-15}$.

```matlab
hold on;
```
**Spiegazione:** Permette di sovrapporre più curve sullo stesso grafico (es. Leja vs Equispaziati).

```matlab
legend;
```
**Spiegazione:** Mostra la legenda per distinguere le curve.

---

## PARTE 3: Interpretazione dei Grafici (Cosa dire all'Esame)

### 1. Grafico Costanti di Lebesgue (Alg 1 vs Alg 2)
<p align="center">
    <img src="LEB.jpg" alt="My figure" width="650"/>
</p>

*   **Visivamente:** Le due curve sono identiche e basse.
*   **Significato:** Entrambi gli algoritmi funzionano correttamente. La costante cresce lentamente (logaritmicamente), il che significa che l'interpolazione di Leja è **stabile** e ben condizionata.

### 2. Grafico Tempo Computazionale
<p align="center">
    <img src="TEMPO.jpg" alt="My figure" width="650"/>
</p>

*   **Visivamente:** Curva Blu (Alg 1) vs Rossa (Alg 2).
*   **Nota Importante:** Potresti vedere l'Alg 1 più veloce dell'Alg 2 per gradi bassi ($d < 50$).
*   **Spiegazione Tecnica:** L'Alg 1 fa solo moltiplicazioni (veloci). L'Alg 2 calcola molti `cosenos` (lenti). Tuttavia, asintoticamente (per gradi altissimi), l'Alg 2 sarebbe preferibile per stabilità numerica.

### 3. Grafico Errore (Leja vs Equispaziati)
<p align="center">
    <img src="ERRORE.jpg" alt="My figure" width="650"/>
</p>

*   **Visivamente:** La curva magenta (Equi) scende e poi risale esplodendo. La verde (Leja) scende sempre fino a zero.
*   **Significato:** Questo grafico dimostra la superiorità dei nodi di Leja. I nodi equispaziati soffrono del **Fenomeno di Runge** a causa della singolarità vicina ($x=1.3$). I nodi di Leja, accumulandosi ai bordi, neutralizzano questo problema e garantiscono la convergenza.
