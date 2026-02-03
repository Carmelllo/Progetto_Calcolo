% MAIN_EXPERIMENT_PROGETTO_1
% Complete script for Progetto 1: Leja Nodes and Lebesgue Constants
clearvars; close all; clc

% 1. Setup parameters
M_leja = 100000;
x_leja = linspace(-1, 1, M_leja)';  
M_eval = 2001;
x_eval = linspace(-1, 1, M_eval)'; 
degrees = 1:50;
f_test = @(x) 1./(x - 1.3);

% Preallocate arrays for results
L_dlp = zeros(length(degrees), 1);
L_dlp2 = zeros(length(degrees), 1);
time_dlp = zeros(length(degrees), 1);
time_dlp2 = zeros(length(degrees), 1);
err_leja = zeros(length(degrees), 1);
err_equi = zeros(length(degrees), 1);

%% Experiments
for ii = 1:length(degrees)
    d = degrees(ii);
    n = d + 1;
    
    tic;
    z1 = dlp(x_leja, d);
    time_dlp(ii) = toc;

    L_dlp(ii) = leb_con(z1', x_eval);
    
    tic;
    z2 = dlp2(x_leja, d);
    time_dlp2(ii) = toc;
    
    L_dlp2(ii) = leb_con(z2', x_eval);
    
    V_leja = zeros(n, n);
    V_eval = zeros(M_eval, n);
    for j = 1:n
        V_leja(:, j) = cos((j-1) * acos(z2));
        V_eval(:, j) = cos((j-1) * acos(x_eval));
    end
    coeff_leja = V_leja \ f_test(z2);
    p_leja = V_eval * coeff_leja;
    err_leja(ii) = max(abs(f_test(x_eval) - p_leja));
    
    zeq = linspace(-1, 1, n)';
    V_eq = zeros(n, n);
    for j = 1:n
        V_eq(:, j) = cos((j-1) * acos(zeq));
    end
    coeff_eq = V_eq \ f_test(zeq);
    p_eq = V_eval * coeff_eq;
    err_equi(ii) = max(abs(f_test(x_eval) - p_eq));
end

%% Visualization

figure('Name', 'Costanti di Lebesgue');
semilogy(degrees, L_dlp, '-b', 'LineWidth', 1.5, 'DisplayName', 'Leja (Alg 1)');
hold on;
semilogy(degrees, L_dlp2, '--r', 'LineWidth', 1.5, 'DisplayName', 'Leja (Alg 2)');
grid on; xlabel('Grado d'); ylabel('\Lambda_n');
title('Confronto delle costanti di Lebesgue'); legend;

figure('Name', 'Tempo computazionale');
plot(degrees, time_dlp, '-b', 'LineWidth', 1.5, 'DisplayName', 'Alg 1 (Produttoria)');
hold on;
plot(degrees, time_dlp2, '-r', 'LineWidth', 1.5, 'DisplayName', 'Alg 2 (LU)');
grid on; xlabel('Grado d'); ylabel('Tempo (s)');
title('Tempo computazionale (Alg 1 vs Alg 2)'); legend;

figure('Name', 'Errore di Interpolazione');
semilogy(degrees, err_leja, '-g', 'LineWidth', 1.5, 'DisplayName', 'Leja (dlp2)');
hold on;
semilogy(degrees, err_equi, '-m', 'LineWidth', 1.5, 'DisplayName', 'Equispaziati');
grid on; xlabel('Grado d'); ylabel('Errore Massimo');
title('Errore di Interpolazione per la funzione: f(x) = 1/(x - 1.3)'); legend;

disp('Esperimenti finiti.');