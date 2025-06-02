function parte_imag = kkimbook(omega, parte_re)
% Calcula la parte imaginaria del índice de refracción/función dieléctrica
% a partir de la parte real n(omega), usando la relación KK.

% Entradas:
% omega - vector de energía en eVs (debe estar equiespaciado)
% nreal - vector de la parte real del índice de refracción/funcion
% dieléctrica

% Salida:
% kimag - vector estimado de la parte imaginaria del índice de
% refracción/función dieléctrica

    % Asegurar que sean vectores fila
    if size(omega,1) > size(omega,2)
        omega = omega';
    end
    if size(parte_re,1) > size(parte_re,2)
        parte_re = parte_re';
    end

    g = length(omega);                  % número de puntos
    parte_imag = zeros(size(parte_re));         % tamaño de la salida
    a = zeros(size(parte_re));             % acumulador izquierdo
    b = zeros(size(parte_re));             % acumulador derecho
    deltaomega = omega(2) - omega(1);   % paso (se asume constante)
    
    % Primer punto (excluye omega(1))
    for k = 2:g
        b(1) = b(1) + (parte_re(k)-1) / (omega(k)^2 - omega(1)^2);
    end
    parte_imag(1) = (-2 * deltaomega * b(1) * omega(1) / pi) ;

    % Último punto (excluye omega(g))
    for k = 1:g-1
        a(g) = a(g) + (parte_re(k)-1) / (omega(k)^2 - omega(g)^2);
    end
    parte_imag(g) = (-2 * deltaomega * a(g) * omega(g)  / pi) ;

    % Puntos intermedios
    for j = 2:g-1

        % Suma desde k = 1 hasta j-1 (antes del punto j)
        for k = 1:j-1
            a(j) = a(j) + (parte_re(k)-1) / (omega(k)^2 - omega(j)^2);
        end

         % Suma desde k = j+1 hasta g (después del punto j)
        for k = j+1:g
            b(j) = b(j) + (parte_re(k)-1) / (omega(k)^2 - omega(j)^2);
        end
        parte_imag(j) = (-2 * deltaomega * (a(j) + b(j)) * omega(j) / pi) ;
    end
end