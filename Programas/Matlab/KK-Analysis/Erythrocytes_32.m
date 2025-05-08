% Define the file path
n_filePath = 'C:\Users\Poh\Documents\Tesis\Refractive index\32_n_ReN_Erythrocytes.csv';

% Read the data into a table, skipping the first 3 rows
n_dataTable = readtable(n_filePath, 'HeaderLines', 1);

% Access the data as a table
%disp(dataTable);
% 
%
lambda = n_dataTable{1:145, 1};         % en nanómetros
n      = n_dataTable{1:145, 2};  % segunda columna
%% 
%Graficar datos iniciales
% scatter(lambda,n)

%%
%Hacer interpolacion de los datos
lambda_interp = linspace(min(lambda), max(lambda), 2000);
n_interp = interp1(lambda, n, lambda_interp, 'spline');

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

%%
%Aplicar KK
k_eri = kkimbook_refractive_index(omega, n_interp);
%%
%Graficar datos obtenidos
% plot(lambda_interp,n_interp)
% hold on
% plot(lambda_interp,n_eri)
% hold on
% legend('k','kOroKK')
%%
[refin, imfin]= selfconsbook_refractive_index(omega, n_interp, k_eri, 500,1);
%%
figure;
subplot(2,1,1)
plot(lambda_interp, n_interp, 'b',lambda_interp, refin, 'r--')
xlabel('\lambda [nm]'); ylabel('n');
legend('n original','n estimado KK'); title('Parte real del índice');

subplot(2,1,2)
plot(lambda_interp, imfin, 'r--')
xlabel('\lambda [nm]'); ylabel('k');
legend('k estimado KK'); title('Parte imaginaria del índice');