function kimag = kkimbook_refractive_index(omega, nreal)
% Calcula la parte imaginaria del índice de refracción k(omega)
% a partir de la parte real n(omega), usando la relación KK.

% Entradas:
% omega en rad/s (debe estar equiespaciado)
% nreal - vector de la parte real del índice de refracción

% Salida:
% kimag - vector estimado de la parte imaginaria del índice de refracción

    % Asegurar que sean vectores fila
    if size(omega,1) > size(omega,2)
        omega = omega';
    end
    if size(nreal,1) > size(nreal,2)
        nreal = nreal';
    end

    g = length(omega);                  % número de puntos
    kimag = zeros(size(nreal));         % tamaño de la salida
    a = zeros(size(nreal));             % acumulador izquierdo
    b = zeros(size(nreal));             % acumulador derecho
    deltaomega = omega(2) - omega(1);   % paso (se asume constante)
    
    % Primer punto (excluye omega(1))
    for k = 2:g
        b(1) = b(1) + (nreal(k)-1) / (omega(k)^2 - omega(1)^2);
    end
    kimag(1) = -2 * deltaomega * b(1) * omega(1) / pi ;

    % Último punto (excluye omega(g))
    for k = 1:g-1
        a(g) = a(g) + (nreal(k)-1) / (omega(k)^2 - omega(g)^2);
    end
    kimag(g) = -2 * deltaomega * a(g) * omega(g) / pi ;

    % Puntos intermedios
    for j = 2:g-1

        % Suma desde k = 1 hasta j-1 (antes del punto j)
        for k = 1:j-1
            a(j) = a(j) + (nreal(k)-1) / (omega(k)^2 - omega(j)^2);
        end

         % Suma desde k = j+1 hasta g (después del punto j)
        for k = j+1:g
            b(j) = b(j) + (nreal(k)-1)/ (omega(k)^2 - omega(j)^2);
        end
        kimag(j) = -2 * deltaomega * (a(j) + b(j)) * omega(j) / pi ;
    end
end