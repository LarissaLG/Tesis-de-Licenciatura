% Define the file path
k_filePath = 'C:\Users\Poh\Documents\Tesis\Refractive index\28.7_k_ImN_Erythrocytes.csv';
n_filePath = 'C:\Users\Poh\Documents\Tesis\Refractive index\28.7_n_ReN_Erythrocytes.csv';

% Read the data into a table, skipping the first 3 rows
k_dataTable = readtable(k_filePath, 'HeaderLines', 1);
n_dataTable = readtable(n_filePath, 'HeaderLines', 1);

% Access the data as a table
%disp(dataTable);
% 
%
n_lambda = n_dataTable{1:90, 1};         % en nanómetros
k_lambda = k_dataTable{1:62, 1};         % en nanómetros
n      = n_dataTable{1:90, 2};  % segunda columna
k      = k_dataTable{1:62, 2}; % segunda columna, distinto rango

%% 
%Graficar datos iniciales
%scatter(n_lambda,n)
% hold on
%scatter(k_lambda,k)

%%
%Hacer interpolacion de los datos
lambda_interp = linspace(min(k_lambda), max(k_lambda), 2000);
n_interp = interp1(n_lambda, n, lambda_interp, 'spline');
k_interp = interp1(k_lambda, k, lambda_interp, 'spline');

%% 
%Graficar datos interpolados
%plot(lambda_interp,n_interp)
% hold on
%plot(lambda_interp,k_interp)

%%
%Convertir lambda a omega
c = 2.9979e8;                          % velocidad de la luz en m/s
lambda_m = lambda_interp * 1e-9;           % convertir nm a m
omega = 2 * pi * c ./ lambda_m;        % calcular omega (rad/s)

%%
%Aplicar KK
k_eri = kkimbook_refractive_index(omega, n_interp);
n_eri = kkrebook_refractive_index(omega, k_interp);

%%
%Graficar datos obtenidos
plot(lambda_interp,n_interp)
hold on
plot(lambda_interp,n_eri)
hold on
legend('k','kOroKK')

%%
%Graficar datos obtenidos
plot(lambda_interp,k_interp)
hold on
plot(lambda_interp,k_eri)
hold on
legend('k','kOroKK')
%%
[refin, imfin]= selfconsbook_refractive_index(omega, n_interp, k_interp, 10, 1);
%%
figure;
subplot(2,1,1)
plot(lambda_interp, n_interp, 'b',lambda_interp, refin, 'r--')
xlabel('\lambda [nm]'); ylabel('n');
legend('n original','n estimado KK'); title('Parte real del índice');

subplot(2,1,2)
plot(lambda_interp, k_interp, 'b',lambda_interp, imfin, 'r--')
xlabel('\lambda [nm]'); ylabel('k');
legend('k original','k estimado KK'); title('Parte imaginaria del índice');