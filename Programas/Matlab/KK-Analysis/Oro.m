% Define the file path
filePath = 'C:\Users\Poh\Documents\GitHub\Servicio-Social\Notebooks Mathematica\Funciones dieléctricas\Oro Johnson.csv';

% Read the data into a table, skipping the first 3 rows
dataTable = readtable(filePath, 'HeaderLines', 1);

% Access the data as a table
%disp(dataTable);
% 
%
lambda_um = dataTable{1:49, 1};         % en micrómetros
lambda_nm = lambda_um * 1000;           % en nanómetros
n      = dataTable{1:49, 2};  % segunda columna
k      = dataTable{51:99, 2}; % segunda columna, distinto rango

%%
%Hacer interpolacion de los datos
lambda_interp = linspace(min(lambda_nm), max(lambda_nm), 10000);
n_interp = interp1(lambda_nm, n, lambda_interp, 'spline');
k_interp = interp1(lambda_nm, k, lambda_interp, 'spline');

%% 
%Graficar datos interpolados
% plot(lambda_interp,n_interp)
% hold on
% plot(lambda_interp,k_interp)

%%
%Convertir lambda a omega
c = 2.9979e8;                          % velocidad de la luz en m/s
lambda_m = lambda_interp * 1e-9;           % convertir nm a m
omega = 2 * pi * c ./ lambda_m;        % calcular omega (rad/s)

% Ordenar omega en orden creciente
[omega_sorted, idx] = sort(omega);
n_sorted = n_interp(idx);

% Aplicar relación KK (usa tu función)
k_KK = kkimbook_refractive_index(omega_sorted, n_sorted);

% Reconvertir a lambda para graficar
lambda_sorted = (2 * pi * c ./ omega_sorted) * 1e6; % en micras

% Graficar
figure;
plot(lambda_interp, k_interp, 'r', 'DisplayName', 'k experimental')
hold on
plot(lambda_interp, k_KK, 'b--', 'DisplayName', 'k calculado por KK')
xlabel('\lambda [\mum]')
ylabel('k')
legend show
title('Comparación de k(\lambda) experimental vs. calculado por Kramers-Kronig')

%%

% Ordenar omega en orden creciente
[omega_sorted, idx] = sort(omega);
k_sorted = k_interp(idx);

% Aplicar relación KK (usa tu función)
n_KK = kkrebook_refractive_index(omega_sorted, k_sorted);

% Reconvertir a lambda para graficar
lambda_sorted = (2 * pi * c ./ omega_sorted) * 1e6; % en micras

% Graficar
figure;
plot(lambda_interp, n_interp, 'r', 'DisplayName', 'k experimental')
hold on
plot(lambda_interp, n_KK, 'b--', 'DisplayName', 'k calculado por KK')
xlabel('\lambda [\mum]')
ylabel('k')
legend show
title('Comparación de k(\lambda) experimental vs. calculado por Kramers-Kronig')


%%
%Aplicar KK
k_oro = kkimbook_refractive_index(omega, n_interp,7);
n_oro = kkrebook_refractive_index(omega, k_interp,7);
%%
%Graficar datos obtenidos
plot(lambda_interp,k_interp)
hold on
plot(lambda_interp,k_oro, 'r--' )
hold on
% plot(lambda_interp,n_oro)
legend('k','kOroKK')

%%
%Graficar datos obtenidos
plot(lambda_interp,n_interp)
hold on
plot(lambda_interp,n_oro)
hold on
% plot(lambda_interp,n_oro)
legend('k','kOroKK')
%%
[refin, imfin]= selfconsbook_refractive_index(omega_sorted, n_sorted, k_KK, 30, 1);
%%
figure;
subplot(2,1,1)
plot(lambda_sorted, n_sorted, 'b',lambda_sorted, refin, 'r--')
xlabel('\lambda [nm]'); ylabel('n');
legend('n original','n estimado KK'); title('Parte real del índice');

subplot(2,1,2)
plot(lambda_sorted, k_sorted, 'b',lambda_sorted, imfin, 'r--')
xlabel('\lambda [nm]'); ylabel('k');
legend('k original','k estimado KK'); title('Parte imaginaria del índice');

%%
% hola2 = kkimbook_refractive_index(lambda_interp, n_interp);
% 
% % Graficar comparación
%     figure;
%     plot(lambda_interp, hola2, 'b', 'LineWidth', 2); hold on;
%     plot(lambda_interp, k_interp, 'c-', 'LineWidth', 2); hold on;
%     plot(lambda_interp, n_interp, 'g', 'LineWidth', 2); hold on;
%     title('Estimación de k(\lambda) mediante Kramers-Kronig');
%     xlabel('\lambda (\mum)');
%     ylabel('k');
%     grid on;
%     legend('k estimado yo (KK)', 'k original');
% 

%%
% Datos experimentales
lambda = [0.1879, 0.1916, 0.1953, 0.1993, 0.2033, 0.2073, 0.2119, 0.2164, ...
          0.2214, 0.2262, 0.2313, 0.2371, 0.2426, 0.249, 0.2551, 0.2616, ...
          0.2689, 0.2761, 0.2844, 0.2924, 0.3009, 0.3107, 0.3204, 0.3315, ...
          0.3425, 0.3542, 0.3679, 0.3815, 0.3974, 0.4133, 0.4305, 0.4509, ...
          0.4714, 0.4959, 0.5209, 0.5486, 0.5821, 0.6168, 0.6595, 0.7045, ...
          0.756, 0.8211, 0.892, 0.984, 1.088, 1.216, 1.393, 1.61, 1.937];

k_exp = [1.188, 1.203, 1.226, 1.251, 1.277, 1.304, 1.35, 1.387, 1.427, ...
         1.46, 1.497, 1.536, 1.577, 1.631, 1.688, 1.749, 1.803, 1.847, ...
         1.869, 1.878, 1.889, 1.893, 1.898, 1.883, 1.871, 1.866, 1.895, ...
         1.933, 1.952, 1.958, 1.948, 1.914, 1.849, 1.833, 2.081, 2.455, ...
         2.863, 3.272, 3.697, 4.103, 4.542, 5.083, 5.663, 6.35, 7.15, ...
         8.145, 9.519, 11.21, 13.78];

n_exp = [1.28, 1.32, 1.34, 1.33, 1.33, 1.3, 1.3, 1.3, 1.3, 1.31, 1.3, ...
         1.32, 1.32, 1.33, 1.33, 1.35, 1.38, 1.43, 1.47, 1.49, 1.53, ...
         1.53, 1.54, 1.48, 1.48, 1.5, 1.48, 1.46, 1.47, 1.46, 1.45, ...
         1.38, 1.31, 1.04, 0.62, 0.43, 0.29, 0.21, 0.14, 0.13, 0.14, ...
         0.16, 0.17, 0.22, 0.27, 0.35, 0.43, 0.56, 0.92];

% Convertir lambda (μm) a frecuencia angular omega (rad/s)
c = 299792458e6; % velocidad de la luz en μm/s
omega = 2 * pi * c ./ lambda; % rad/s

% Ordenar omega de menor a mayor (necesario para integración)
[omega_sorted, idx] = sort(omega);
k_sorted = k_exp(idx);
n_sorted = n_exp(idx);

% Parámetro alpha
alpha = 4;

% Aplicar Kramers-Kronig inversa: obtener n desde k
g = length(omega_sorted);
reN = zeros(1, g);
deltaomega = omega_sorted(2) - omega_sorted(1);

for j = 1:g
    suma = 0;
    for k = [1:j-1, j+1:g]
        suma = suma + k_sorted(k)*(omega_sorted(k)^(2*alpha+1)) / ...
                     (omega_sorted(k)^2 - omega_sorted(j)^2);
    end
    reN(j) = (2/pi * deltaomega * suma * omega_sorted(j)^(-2 *alpha)) + 1;
end

% Aplicar Kramers-Kronig inversa: obtener k desde n
g = length(omega_sorted);
imN = zeros(1, g);
deltaomega = omega_sorted(2) - omega_sorted(1);

for j = 1:g
    suma = 0;
    for k = [1:j-1, j+1:g]
        suma = suma + (n_sorted(k)-1)*(omega_sorted(k)^(2*alpha)) / ...
                     (omega_sorted(k)^2 - omega_sorted(j)^2);
    end
    imN(j) = -(2/pi * deltaomega * suma * omega_sorted(j)^(1-2 *alpha));
end

% Graficar en lambda (longitud de onda)
lambda_sorted = lambda(idx);

figure;
plot(lambda_sorted, n_sorted, 'b', 'LineWidth', 1.5); hold on;
plot(lambda_sorted, reN, 'r--', 'LineWidth', 1.5);
% plot(lambda_sorted, k_sorted, 'g', 'LineWidth', 1.5); hold on;
plot(lambda_sorted, imN, 'm--', 'LineWidth', 1.5);
xlabel('Longitud de onda (μm)');
ylabel('Índice de refracción real');
legend('n experimental', 'n desde k (KK)','k experimental','k desde n (KK)');
title('Kramers-Kronig: cálculo de n desde k');
grid on;
