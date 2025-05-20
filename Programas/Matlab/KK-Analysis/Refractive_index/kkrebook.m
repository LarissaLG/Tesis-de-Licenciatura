function nreal = kkrebook(omega, kimag)
% Calcula la parte real del índice de refracción/función dieléctrica a partir
% de su parte imaginaria mediante las relaciones de Kramers-Kronig.
% Entradas:
%   omega  - vector de energia en eVs (debe estar equiespaciado)
%   kimag  - vector de la parte imaginaria del índice de refracción/funcion
% dieléctrica
%
% Salida:
%   nreal  - vector estimado de la parte real del índice de refracción/funcion
% dieléctrica

    % Asegura que omega e imN sean vectores fila
    if size(omega,1) > size(omega,2)
        omega = omega';
    end
    if size(kimag,1) > size(kimag,2)
        kimag = kimag';
    end

    g = length(omega);                    % número de puntos
    nreal = zeros(1, g);                  % inicializa salida
    a = zeros(1, g);                      % acumulador izquierdo
    b = zeros(1, g);                      % acumulador derecho
    deltaomega = omega(2) - omega(1);     % paso (se asume constante)

    % Primer punto (excluye omega(1))
    for k = 2:g
        b(1) = b(1)+ kimag(k)* omega(k)/ (omega(k)^2 - omega(1)^2);
    end
    nreal(1) = (2/pi * deltaomega * b(1))+1;

    % Último punto (excluye omega(g))
    for k = 1:g-1
        a(g) = a(g)+ kimag(k)*omega(k)/ (omega(k)^2 - omega(g)^2);
    end
    nreal(g) = (2/pi * deltaomega * a(g))+1;

    % Puntos intermedios
    for j = 2:g-1
  
        % Suma desde k = 1 hasta j-1 (antes del punto j)
        for k = 1:j-1
            a(j) = a(j) + kimag(k)* omega(k)/ (omega(k)^2 - omega(j)^2);
        end

        % Suma desde k = j+1 hasta g (después del punto j)
        for k = j+1:g
            b(j) = b(j) + kimag(k)*omega(k)/ (omega(k)^2 - omega(j)^2);
        end

        % Parte real aproximada en omega(j)
        nreal(j) = (2/pi * deltaomega * (a(j) + b(j)))+1;
    end
end
