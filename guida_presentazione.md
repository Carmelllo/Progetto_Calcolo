# Guida Approfondita per la Presentazione (20 Minuti)

Questo documento è strutturato per supportare una spiegazione orale dettagliata e tecnica. Ogni sezione include la **Logica Matematica**, la **Sintassi MATLAB**, e i **Punti Chiave per l'Esposizione**.

---

## 1. Introduzione Teorica (Concetti Fondamentali)

**Da dire all'inizio:**
> "Il problema fondamentale dell'interpolazione polinomiale su nodi equispaziati è l'instabilità numerica, nota come **Fenomeno di Runge**. All'aumentare del grado $N$, l'errore ai bordi esplode esponenzialmente.
> L'obiettivo di questo progetto è implementare e analizzare i **Punti di Leja Approssimati**, una distribuzione di nodi che garantisce stabilità (crescita logaritmica della costante di Lebesgue) e che può essere calcolata in modo efficiente sfruttando la proprietà di essere una successione annidata (cosa che i nodi di Chebyshev non sono)."

---

## 2. Analisi Dettagliata dei File `.m`

### File 1: `dlp.m` (Algoritmo Greedy / Iterativo)

**Concetto:** Questo algoritmo implementa la definizione formale dei punti di Leja. È un approccio "Greedy" (ingordo): cerca l'ottimo locale ad ogni passo.

```matlab
function dlp = dlp(x, d)
    n = d + 1;
    dlp = zeros(n, 1);  % Preallocazione per Efficienza
    dlp(1) = x(1);      % Inizializzazione Arbitraria (-1)
    
    for k = 2:n
        prod_val = ones(length(x), 1);
        
        % Calcolo produttoria distanze dai nodi precedenti
        for j = 1:k-1
            prod_val = prod_val .* abs(x - dlp(j));
        end
        
        % Scelta del punto che massimizza il prodotto
        [~, idx] = max(prod_val);
        dlp(k) = x(idx);
    end
end
```

**Spiegazione Approfondita per l'Orale:**
1.  **Input:** La `x` non è l'intervallo continuo, ma una "mesh discretizzata" molto fitta (es. 100.000 punti). Stiamo scegliendo i candidati migliori da questo insieme finito.
2.  **Il Cuore Matematico (`prod_val`):**
    *   Ad ogni passo $k$, costruiamo un "potenziale" su tutta la mesh.
    *   La formula `prod_val = prod_val .* abs(x - dlp(j))` calcola $\prod_{j=0}^{k-1} |x - \xi_j|$.
    *   Fisicamente, questo è analogo a trovare dove il potenziale elettrico è minimo (o dove la repulsione è minima) tra cariche dello stesso segno. I punti cercano di stare il più lontani possibile.
3.  **Complessità:** Questo algoritmo ha due cicli annidati (`for k` e `for j`) dentro operazioni su vettori lunghi $M$. La complessità è circa $O(M \cdot N^2)$. È lento per gradi molto alti ($N > 100$), motivo per cui esiste `dlp2`.
4.  **Sintassi MATLAB:**
    *   `.*` (Punto-Per): Sottolinea che questa è una moltiplicazione elemento-per-elemento, non matriciale.
    *   `max(Val)`: Restituisce sia il massimo che l'indice. L'indice è ciò che ci serve per localizzare il punto nella mesh.

---

### File 2: `dlp2.m` (Algoritmo Matriciale / Algebra Lineare)

**Concetto:** Sfrutta la relazione tra i punti di Leja e il determinante della matrice di Vandermonde.

```matlab
function dlp = dlp2(x, d)
    n = d + 1;
    M = length(x);
    V = zeros(M, n);
    
    % Costruzione Matrice di Vandermonde-Chebyshev
    for j = 1:n
        V(:, j) = cos((j-1) * acos(x));
    end
    
    % Fattorizzazione LU con Pivoting
    [~, ~, p] = lu(V, 'vector');
    
    dlp = x(p(1:n));
end
```

**Spiegazione Approfondita per l'Orale:**
1.  **Cambio di Paradigma:** Invece di misurare distanze geometriche, cerchiamo le righe della matrice che sono "più linearmente indipendenti".
2.  **La Base di Chebyshev (`cos(k*acos(x))`):**
    *   Non usiamo la base monomiale $1, x, x^2$ perché la matrice di Vandermonde classica è malcondizionata (gli elementi diventano piccolissimi o grandissimi).
    *   I polinomi di Chebyshev oscillano sempre tra -1 e 1. Questo garantisce stabilità numerica nel calcolo della fattorizzazione.
3.  **Il Ruolo di LU:**
    *   La Gaussian Elimination con **Partial Pivoting** scambia le righe per portare l'elemento più grande ("pivot") sulla diagonale.
    *   Matematicamente, scegliere il pivot più grande equivale a massimizzare il volume del parallelepipedo formato dai vettori riga.
    *   Le righe scelte per prime corrispondono ai punti che massimizzano il determinante, che è equivalente alla condizione di Leja.
4.  **Efficienza:** Questo metodo delega il lavoro pesante alle librerie LAPACK (ottimizzate in Fortran/C) usate da MATLAB per `lu`. È molto più veloce per gradi elevati.

---

### File 3: `leb_con.m` (Indicatori di Qualità)

**Concetto:** Come misuriamo se i punti scelti sono "buoni"? Usiamo la Costante di Lebesgue.

```matlab
function L = leb_con(z, x)
    % ... inizializzazione ...
    for i = 1:n+1
        % Calcolo i-esimo polinomio base di Lagrange
        li = ones(m, 1);
        for j = [1:i-1, i+1:n+1] % Salta l'indice i
            li = li .* (x - z(j)) / (z(i) - z(j));
        end
        lambda = lambda + abs(li);
    end
    L = max(lambda);
end
```

**Spiegazione Approfondita per l'Orale:**
1.  **Significato di $\Lambda_n(x)$:** La funzione di Lebesgue $\Lambda_n(x) = \sum |l_i(x)|$ rappresenta il fattore di amplificazione dell'errore.
    *   Se misuriamo i dati $y$ con un piccolo errore $\epsilon$, l'errore sul polinomio interpolante sarà limitato da $L \cdot \epsilon$.
2.  **Perché `max(lambda)`?** Vogliamo proteggerci dal "caso peggiore". La costante $L$ è il picco massimo di amplificazione sull'intervallo.
3.  **Confronto Teorico:**
    *   Equispaziati: $L \sim 2^N$. Disastroso.
    *   Leja/Chebyshev: $L \sim \frac{2}{\pi} \log(N)$. Ottimale (crescita lentissima).

---

### File 4: `main_experiment.m` (Analisi Risultati)

**Concetto:** Lo script che orchestra l'esperimento, confronta i metodi e visualizza i dati.

**Punti Critici del Codice:**

1.  **`clearvars; close all;`**:
    > "Iniziare sempre con un ambiente pulito garantisce la riproducibilità. Evita che variabili 'fantasma' di sessioni precedenti falsino i risultati."

2.  **`f_test = @(x) 1./(x - 1.3);`**:
    > "Abbiamo scelto una funzione di Runge con una singolarità in $x=1.3$. Essendo vicina al bordo dell'intervallo ($x=1$), questa funzione mette a dura prova la stabilità dell'interpolazione ai bordi."

3.  **Risoluzione del Sistema (`V \ y`):**
    ```matlab
    coeff_leja = V_leja \ f_test(z2);
    ```
    > **Importante:** "Non calcoliamo l'inversa $V^{-1}$. In analisi numerica, l'inversa è costosa e instabile. L'operatore backslash usa la fattorizzazione (Gauss/QR) per risolvere il sistema direttamente. È lo standard industriale."

4.  **Visualizzazione (`semilogy`):**
    > "Usiamo la scala logaritmica per l'asse Y perché stiamo confrontando errori che variano di 15 ordini di grandezza (da $10^0$ a $10^{-15}$). In scala lineare, l'errore di Leja sembrerebbe zero e non apprezzeremmo la sua precisione."

---

## 3. Interpretazione dei Grafici (Cosa dire mostrando le figure)

### Grafico 1: Tempi di Calcolo
*   **Osservazione:** "Notiamo che per gradi bassi ($N < 50$), l'algoritmo Greedy (`dlp`) è competitivo. Tuttavia, teoricamente, l'algoritmo matriciale (`dlp2`) scala meglio asintoticamente."
*   **Dettaglio:** "La curva `dlp2` potrebbe apparire più piatta perché gran parte del tempo è spesa nell'allocazione della matrice, mentre `dlp` cresce quadraticamente."

### Grafico 2: Costanti di Lebesgue
*   **Osservazione:** "Le curve per `dlp` e `dlp2` sono sovrapposte e crescono molto lentamente (logaritmicamente). Questo conferma che entrambi gli algoritmi estraggono punti con la stessa qualità asintotica dei nodi di Chebyshev."
*   **Confronto:** "Se avessimo plottato i nodi equispaziati, la curva sarebbe schizzata verticalmente fuori dal grafico già per $N=20$."

### Grafico 3: Errore di Interpolazione
*   **Analisi:** "L'errore sui nodi di Leja decresce fino alla precisione di macchina ($\approx 10^{-14}$), dimostrando la convergenza spettrale per funzioni analitiche."
*   **Equispaziati:** "Al contrario, con nodi equispaziati, l'errore esplode ai bordi a causa del fenomeno di Runge, rendendo l'interpolazione inutile per gradi alti."

---

## 4. Domande Probabili e Risposte

**Q: Perché usiamo una "mesh fitta" invece dell'intervallo continuo?**
A: I punti di Leja reali sono definiti sul continuo nel piano complesso (insiemi di capacità). Discretizzare l'intervallo in una mesh $K$ permette di trasformare un problema di ottimizzazione continua in uno di ricerca combinatoria, molto più facile da risolvere al calcolatore.

**Q: Perché la matrice di Vandermonde è tipicamente malcondizionata?**
A: Le colonne $1, x, x^2, \dots$ diventano quasi parallele (linearmente dipendenti) per $x \in [0,1]$ quando il grado è alto. Usare la base ortogonale di Chebyshev ("Vandermonde generalizzata") mantiene le colonne ben separate angolarmente, migliorando drasticamente il condizionamento.

**Q: Qual è il vantaggio pratico dei punti di Leja rispetto a Chebyshev?**
A: I punti di Chebyshev sono definiti rigidamente per ogni grado $N$. Se decido di passare da grado 10 a grado 11, devo ricalcolare *tutti* i punti e *tutti* i valori della funzione. La sequenza di Leja è **annidata**: per passare da $N$ a $N+1$ aggiungo solo un nuovo punto, riutilizzando tutto il lavoro precedente. Questo è fondamentale per algoritmi adattivi (es. "continua ad aggiungere punti finché l'errore non scende sotto la soglia").
