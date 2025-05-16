function reN = sskkrebook_refractive_index(omega, imN, reN1, omega1)
% sskkrebook - Calcula la parte real del índice de refracción a partir
% de su parte imaginaria mediante las relaciones de Kramers-Kronig
% sustractivas
%
% Entradas:
%   omega  - vector de frecuencias (debe estar equiespaciado)
%   imN  - vector de la parte imaginaria del índice de refracción
%
% Salida:
%   reN  - vector estimado de la parte real del índice de refracción

    % Asegura que omega e imN sean vectores fila
    if size(omega,1) > size(omega,2)
        omega = omega';
    end
    if size(imN,1) > size(imN,2)
        imN = imN';
    end

    g = length(omega);                     % número de puntos
    reN = zeros(1, g);                  % inicializa salida
    a = zeros(1, g);                      % acumulador izquierdo
    b = zeros(1, g);                      % acumulador derecho
    deltaomega = omega(2) - omega(1);     % paso de frecuencia (se asume constante)

    % Encuentra el índice del omega más cercano a omega1
    [~, idx] = min(abs(omega - omega1));
    
    % % Pon el valor de reN1
    % reN(idx) = reN1;


    % ---- Primer punto: j = 1 ----
    for k = 2:g
        if abs(omega(k) - omega1) < 1e-3
            continue; % evita división por cero si omega(k) es omega1
        end
        b(1) = b(1)+ imN(k)* omega(k)/ ((omega(k)^2 - omega(1)^2) * (omega(k)^2 - omega1^2));
    end
    reN(1) = (2/pi * deltaomega * b(1) * (omega(1)^2 - omega1^2)) + reN1;

    % ---- Último punto: j = g ----
    for k = 1:g-1
        if abs(omega(k) - omega1) < 1e-3
            continue; % evita división por cero si omega(k) es omega1
        end
        a(g) = a(g)+ imN(k)*omega(k)/ ((omega(k)^2 - omega(g)^2) * (omega(k)^2 - omega1^2));
    end
    reN(g) = (2/pi * deltaomega * a(g) * (omega(g)^2 - omega1^2) ) + reN1;

    % ---- Puntos intermedios: j = 2 hasta g-1 ----
    for j = 2:g-1
        
        if abs(omega(k) - omega1) < 1e-3
            continue; % evita división por cero si omega(k) es omega1
        end 

        % Suma desde k = 1 hasta j-1 (antes del punto j)
        for k = 1:j-1

            if abs(omega(k) - omega1) < 1e-3
                continue; % evita división por cero si omega(k) es omega1
            end 

            a(j) = a(j) + imN(k)* omega(k)/ ((omega(k)^2 - omega(j)^2)  * (omega(k)^2 - omega1^2));
        end

        % Suma desde k = j+1 hasta g (después del punto j)
        for k = j+1:g

            if abs(omega(k) - omega1) < 1e-5
                continue; % evita división por cero si omega(k) es omega1
            end 

            b(j) = b(j) + imN(k)*omega(k) / ((omega(k)^2 - omega(j)^2)  * (omega(k)^2 - omega1^2));
        end

        % Parte real aproximada en omega(j)
        reN(j) = (2/pi * deltaomega * (a(j) + b(j)) * (omega(j)^2 - omega1^2)) + reN1;
    end
end
