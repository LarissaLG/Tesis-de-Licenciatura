function [epsilon, n_complex] = drude_model(omega, omega_p, gamma)
% DRUDE_MODEL - Calcula la epsilon y el índice de refracción usando el modelo de Drude
%
% Entradas:
%   omega    -  (rad/s)
%   omega_p  - Frecuencia de plasma (rad/s)
%   gamma    - Damping  (rad/s)
%
% Salidas:
%   epsilon     - Función dieléctrica(omega)
%   n_complex   - Índice de refracción: n + i*k

    % Modelo de drude
    epsilon = 1 - (omega_p^2) ./ (omega.^2 + 1i * gamma .* omega);

    % Complex refractive index
    n_complex = sqrt(epsilon);
end

%%
% Constantes
c = 3e8; 
hbar = 6.582119569 * 1e-16;

% % Rango de longitud de onda
% lambda_um = linspace(0.2, 1, 1000); % micrómetros
% lambda_nm = lambda_um*1000; % nanómetros 
% lambda_m = lambda_um * 1e-6; %metros

% Convert to angular frequency (rad/s)
% omega = 2 * pi *c ./lambda_m;
omega = linspace(0.5, 8, 1000);

% Drude parameters (example values for aluminim)
omega_p = 4.3; %13.142/hbar; % rad/s
gamma =  0.15;%0.197/hbar;   % rad/s

[eps, n_complex] = drude_model(omega, omega_p, gamma);
n = real(eps);%real(n_complex);
k = imag(eps);%imag(n_complex);

% Plot
figure;
plot(omega, n, 'b', 'LineWidth', 2); hold on;
plot(omega, k, 'r', 'LineWidth', 2);
xlabel('\lambda (nm)');
ylabel('n, k');
legend('n(\lambda)', 'k(\lambda)');
title('Drude Model: Complex Refractive Index');
grid on;

%%
% % Ordenar omega de menor a mayor
% [omega_sorted, idx] = sort(omega);
% 
% % Reordenar k de acuerdo al nuevo orden de omega
% k_sorted = k(idx);
% n_sorted = n(idx);
% 
% %Para regresar a lambda
% lambda_sorted = 2 * pi * 3e8 ./ omega_sorted; % en metros
% lambda_nm_sorted = lambda_sorted * 1e9;


%%
%Aplicar KK
k_oro = kkimbook_refractive_index(omega, n);
n_oro = kkrebook_refractive_index(omega, k);

%%
% [lambda_nm_sorted_up, idx2] = sort(lambda_nm_sorted);
% n_oro_sorted = n_oro(idx2);
% k_oro_sorted = k_oro(idx2);

%%
plot(omega,k_oro, 'r--')
hold on
plot(omega,n_oro,'b--')
hold on
plot(omega, n, 'b', 'LineWidth', 2); hold on;
plot(omega, k, 'r', 'LineWidth', 2);
%%
%Graficar datos obtenidos
plot(omega_sorted,k_sorted, 'r')
hold on
plot(omega_sorted,k_oro_sorted,'r--')
hold on
plot(omega_sorted,n_sorted, 'b')
hold on
plot(omega_sorted,n_oro_sorted,'b--')
legend('n','nKK','k','kKK')

%%
[refin, imfin]= selfconsbook_refractive_index(omega, n_oro, k_oro, 200,0.5);
% refin_sorted = refin(idx2);
% imfin_sorted = imfin(idx2);

%%
%Graficar datos obtenidos
%%
plot(omega,k_oro, 'r--')
hold on
plot(omega,n_oro,'b--')
hold on
plot(omega,imfin, 'r:','LineWidth', 2)
hold on
plot(omega,refin,'b:','LineWidth', 2)
hold on
plot(omega, n, 'b'); hold on;
plot(omega, k, 'r');
%%
figure;
subplot(2,1,1)
plot(lambda_nm, n, 'b',lambda_nm, refin, 'r--')
xlabel('\lambda [nm]'); ylabel('n');
xlim([0 1000])
legend('n original','n estimado KK'); title('Parte real del índice');

subplot(2,1,2)
plot(lambda_nm, k, 'b',lambda_nm, imfin_sorted, 'r--')
xlabel('\lambda [nm]'); ylabel('k');
xlim([0 1000])
legend('k original','k estimado KK'); title('Parte imaginaria del índice');


%% 
% Subtractive K-K Relations

ren= sskkrebook_refractive_index(omega_sorted, k_sorted, 0.8071, 8.7585*1.0e+15);

ren_sorted = ren(idx2);

%%
%Graficar datos obtenidos
% plot(lambda_nm,k, 'r')
% hold on
% plot(lambda_nm,k_oro_sorted,'r--')
% hold on
plot(omega_sorted,n_sorted, 'b')
hold on
plot(omega_sorted,ren,'b--')
legend('n','nKK','k','kKK')
