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
lambda_interp = linspace(min(lambda_nm), max(lambda_nm), 5000);
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

%% 
% Definimos la susceptibilidad eléctrica
% chiRe = n_interp.^ 2 - k_interp.^2 - 1;
% chiIm = 2*n_interp.* k_interp;

%%
%Aplicar KK
k_oro = kkimbook_refractive_index(omega, n_interp);
n_oro = kkrebook_refractive_index(omega, k_interp);
%%
%Graficar datos obtenidos
plot(lambda_interp,k_interp)
hold on
plot(lambda_interp,k_oro)
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
[refin, imfin]= selfconsbook_refractive_index(omega, n_interp, k_interp, 100,1);
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
