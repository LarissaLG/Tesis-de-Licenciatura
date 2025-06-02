function imchi = kkimbookchi(omega, rechi, alpha)
% kkimbook - Calcula la parte imaginaria de la susceptibilidad
% a partir de la parte real mediante las relaciones de Kramers-Kronig.
%
% Entradas:
%   omega  - vector de frecuencias (debe ser equiespaciado)
%   rechi  - vector de la parte real de la susceptibilidad
%   alpha  - parámetro del momento generalizado
%
% Salida:
%   imchi  - vector de la parte imaginaria estimada

    % Asegurar que omega y rechi sean vectores fila
    if size(omega,1) > size(omega,2)
        omega = omega';
    end
    if size(rechi,1) > size(rechi,2)
        rechi = rechi';
    end

    g = length(omega);                     % número de puntos
    imchi = zeros(1, g);                  % inicializa salida
    a = zeros(1, g);                      % acumulador lado izquierdo
    b = zeros(1, g);                      % acumulador lado derecho
    deltaomega = omega(2) - omega(1);    % paso en omega (se supone equiespaciado)

    % ---- Primer punto: j = 1 ----
    beta1 = 0;
    for k = 2:g
        b(1) = beta1 + rechi(k) * omega(k)^(2*alpha) / (omega(k)^2 - omega(1)^2);
        beta1 = b(1);
    end
    imchi(1) = -2/pi * deltaomega * b(1) * omega(1)^(1 - 2*alpha);

    % ---- Último punto: j = g ----
    alpha1 = 0;
    for k = 1:g-1
        a(g) = alpha1 + rechi(k) * omega(k)^(2*alpha) / (omega(k)^2 - omega(g)^2);
        alpha1 = a(g);
    end
    imchi(g) = -2/pi * deltaomega * a(g) * omega(g)^(1 - 2*alpha);

    % ---- Puntos intermedios: j = 2 hasta g-1 ----
    for j = 2:g-1
        alpha1 = 0;
        beta1 = 0;

        for k = 1:j-1
            a(j) = alpha1 + rechi(k) * omega(k)^(2*alpha) / (omega(k)^2 - omega(j)^2);
            alpha1 = a(j);
        end

        for k = j+1:g
            b(j) = beta1 + rechi(k) * omega(k)^(2*alpha) / (omega(k)^2 - omega(j)^2);
            beta1 = b(j);
        end

        imchi(j) = -2/pi * deltaomega * (a(j) + b(j)) * omega(j)^(1 - 2*alpha);
    end
end

