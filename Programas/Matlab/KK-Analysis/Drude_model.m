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
omega = linspace(0.01, 5, 1000);

% Parámetros de Drude 
omega_p = 13.142; %eV
gamma =  0.197; %eV

[eps, n_complex] = drude_model(omega, omega_p, gamma);
n = real(n_complex);
k = imag(n_complex);

epsreal = real(eps);
epsimag = imag(eps);

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
errork = abs(k - kKK);
errorn = abs(n - nKK);

plot(omega, errorn, 'b--')
hold on
plot(omega, errork,'r--')
hold on
xlabel('$\hbar\omega$ (eV)','Interpreter','latex');
legend('Error parte real', 'Error parte imaginaria');


%%
[refin, imfin]= selfconsbook(omega, nKK, k, 1,0.5);
[refin1, imfin1]= selfconsbook(omega, nKK, k, 3,0.5);
%%
%Graficar datos obtenidos
plot(omega,imfin, 'r:','LineWidth', 2)
hold on
plot(omega,refin,'b:','LineWidth', 2)
hold on
plot(omega,imfin1, 'r--')
hold on
plot(omega,refin1,'b--')
hold on
plot(omega, n, 'b'); hold on;
plot(omega, k, 'r');


%%
filename = 'C:\Users\Poh\Documents\GitHub\Tesis-de-Licenciatura\Programas\Matlab\KK-Analysis\selconsbook.gif';
figure; 

frameCount = 1;

for N = 1:1:20
    % Reemplaza esta línea por tu función real
    [refin, imfin] = selfconsbook(omega, nKK, kKK, N, 0.5);

    % Gráfico
    plot(omega, refin, 'b--', omega, imfin, 'r--', omega, n, 'b',omega, k, 'r');
    xlabel('$\hbar\omega$ (eV)','Interpreter','latex');
    title(['N = ', num2str(N)]);
    legend('Parte real KK', 'Parte imaginaria KK','Parte real', 'Parte imaginaria');
    xlim([min(omega), max(omega)]);
    drawnow;

    % Capturar y guardar frame
    frame = getframe(gcf);
    img = frame2im(frame);
    [imind, cm] = rgb2ind(img, 256);

    if frameCount == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.5);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
    end
    frameCount = frameCount + 1;
end





%%
omega1 = 6.0;                     % por ejemplo
nreal1 = interp1(omega, n, omega1);  % del modelo de Drude

kSSKK = sskkimbook(omega,n,omega1,nreal1);
%%
%Graficar datos obtenidos
% plot(omega,kKK, 'r--')
% hold on
% plot(omega,nKK,'b--')
% hold on
plot(omega,kSSKK, 'r:','LineWidth', 2)
hold on
% plot(omega, n, 'b'); hold on;
plot(omega, k, 'r');
legend('kKK', 'nKK','k','n');

%%

nSSKK = sskkrebook2(omega,k,omega1,nreal1);
kKK = kkimbook(omega, nSSKK);


%%
%Graficar datos obtenidos
plot(omega,kKK, 'r--')
hold on
% plot(omega,nKK,'b--')
% hold on
plot(omega,nSSKK, 'b:','LineWidth', 2)
hold on
plot(omega, n, 'b'); hold on;
plot(omega, k, 'r');