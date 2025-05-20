function [epsilon, n_complex] = drude_model(omega, omega_p, gamma)
% DRUDE_MODEL - Calcula la epsilon y el índice de refracción usando el modelo de Drude
%
% Entradas:
%   omega    -  (eV)
%   omega_p  - Frecuencia de plasma (eV)
%   gamma    - Damping  (eV)
%
% Salidas:
%   epsilon     - Función dieléctrica
%   n_complex   - Índice de refracción: n + i*k

    % Modelo de drude
    epsilon = 1 - (omega_p^2) ./ (omega.^2 + 1i * gamma .* omega);

    % Complex refractive index
    n_complex = sqrt(epsilon);
end

%%
% Constantes
hbar = 6.582119569 * 1e-16; % (eV s)

% Vector de energías
omega = linspace(1, 5, 1000);

% Parámetros de Drude 
omega_p = 13.142; %eV
gamma =  0.197; %eV

[eps, n_complex] = drude_model(omega, omega_p, gamma);
n = real(n_complex);
k = imag(n_complex);

epsreal = real(eps);%real(n_complex);
epsimag = imag(eps);%imag(n_complex);

% Plot
figure;
plot(omega, n, 'b', 'LineWidth', 2); hold on;
plot(omega, k, 'r', 'LineWidth', 2);
xlabel('$\hbar\omega$ (eV)','Interpreter','latex');
ylabel('Re(\epsilon), Im(\epsilon)');
legend('Re(\epsilon)', 'Im(\epsilon)');
title('Drude Model: Complex Refractive Index');
grid on;

%%
%Aplicar KK
kKK = kkimbook(omega, n);
nKK = kkrebook(omega, k);

%%
plot(omega, kKK, 'r--')
hold on
plot(omega, nKK,'b--')
hold on
plot(omega, n, 'b', 'LineWidth', 2); hold on;
plot(omega, k, 'r', 'LineWidth', 2);
xlabel('$\hbar\omega$ (eV)','Interpreter','latex');
legend('kKK', 'nKK','n','k');
%%
[refin, imfin]= selfconsbook(omega, nKK, kKK, 200,0.5);
%%
%Graficar datos obtenidos
plot(omega,kKK, 'r--')
hold on
plot(omega,nKK,'b--')
hold on
plot(omega,imfin, 'r:','LineWidth', 2)
hold on
plot(omega,refin,'b:','LineWidth', 2)
hold on
plot(omega, n, 'b'); hold on;
plot(omega, k, 'r');
%%
omega1 = 2.0;                     % por ejemplo
kimag1 = interp1(omega, k, omega1);  % del modelo de Drude

kSSKK = sskkimbook(omega,n,omega1,kimag1);
%%
%Graficar datos obtenidos
plot(omega,kKK, 'r--')
hold on
plot(omega,nKK,'b--')
hold on
plot(omega,kSSKK, 'r:','LineWidth', 2)
hold on
plot(omega, n, 'b'); hold on;
plot(omega, k, 'r');
