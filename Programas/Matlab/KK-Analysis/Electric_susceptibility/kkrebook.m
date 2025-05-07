function rechi = kkrebook(omega, imchi, alpha)
% kkrebook - Calcula la parte real de la susceptibilidad a partir
% de su parte imaginaria mediante las relaciones de Kramers-Kronig.
%
% Entradas:
%   omega  - vector de frecuencias (debe estar equiespaciado)
%   imchi  - vector de la parte imaginaria de la susceptibilidad
%   alpha  - momento generalizado de la transformación
%
% Salida:
%   rechi  - vector estimado de la parte real de la susceptibilidad

    % Asegura que omega e imchi sean vectores fila
    if size(omega,1) > size(omega,2)
        omega = omega';
    end
    if size(imchi,1) > size(imchi,2)
        imchi = imchi';
    end

    g = length(omega);                     % número de puntos
    rechi = zeros(1, g);                  % inicializa salida
    a = zeros(1, g);                      % acumulador izquierdo
    b = zeros(1, g);                      % acumulador derecho
    deltaomega = omega(2) - omega(1);     % paso de frecuencia (se asume constante)

    % ---- Primer punto: j = 1 ----
    beta1 = 0;
    for k = 2:g
        b(1) = beta1 + imchi(k) * omega(k)^(2*alpha + 1) / (omega(k)^2 - omega(1)^2);
        beta1 = b(1);
    end
    rechi(1) = 2/pi * deltaomega * b(1) * omega(1)^(-2*alpha);

    % ---- Último punto: j = g ----
    alpha1 = 0;
    for k = 1:g-1
        a(g) = alpha1 + imchi(k) * omega(k)^(2*alpha + 1) / (omega(k)^2 - omega(g)^2);
        alpha1 = a(g);
    end
    rechi(g) = 2/pi * deltaomega * a(g) * omega(g)^(-2*alpha);

    % ---- Puntos intermedios: j = 2 hasta g-1 ----
    for j = 2:g-1
        alpha1 = 0;
        beta1 = 0;

        % Suma desde k = 1 hasta j-1 (antes del punto j)
        for k = 1:j-1
            a(j) = alpha1 + imchi(k) * omega(k)^(2*alpha + 1) / (omega(k)^2 - omega(j)^2);
            alpha1 = a(j);
        end

        % Suma desde k = j+1 hasta g (después del punto j)
        for k = j+1:g
            b(j) = beta1 + imchi(k) * omega(k)^(2*alpha + 1) / (omega(k)^2 - omega(j)^2);
            beta1 = b(j);
        end

        % Parte real aproximada en omega(j)
        rechi(j) = 2/pi * deltaomega * (a(j) + b(j)) * omega(j)^(-2*alpha);
    end
end
 
%% 
function rechi = drude_rechi(omega, wp, gamma)
% DRUDE_RECHI - Real part of the susceptibility using the Drude model
%
% Inputs:
%   omega - Frequency vector [rad/s]
%   wp    - Plasma frequency [rad/s]
%   gamma - Damping constant [rad/s]
%
% Output:
%   rechi - Real part of the susceptibility χ'(ω)

% Drude model: χ'(ω) = - (wp^2 * (omega.^2 - gamma^2)) / ((omega.^2 + gamma^2).^2)
rechi = - (wp^2 .* (omega.^2 - gamma^2)) ./ ((omega.^2 + gamma^2).^2);

end
%% 

omega = linspace(1e13, 1e16, 1000); % frecuencia en rad/s
wp = 1e16;    % frecuencia de plasma
gamma = 1e14; % tasa de amortiguamiento

rechi = drude_rechi(omega, wp, gamma);
plot(omega, rechi);
xlabel('\omega [rad/s]');
ylabel('Re[\chi(\omega)]');
title('Susceptibilidad real - Modelo de Drude');
grid on;
