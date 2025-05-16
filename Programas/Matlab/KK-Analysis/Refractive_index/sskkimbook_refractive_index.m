function kimag = sskkimbook_refractive_index(omega, nreal, kimag1,omega1)
% Calcula la parte imaginaria del índice de refracción k(omega)
% a partir de la parte real n(omega), usando la relación SSKK.
% omega en rad/s, nreal sin unidad (adimensional)
% kimag es la parte imaginaria del índice de refracción
% k1

    % Asegurar que sean vectores fila
    if size(omega,1) > size(omega,2)
        omega = omega';
    end
    if size(nreal,1) > size(nreal,2)
        nreal = nreal';
    end

    g = length(omega);
    kimag = zeros(size(nreal));
    a = zeros(size(nreal));
    b = zeros(size(nreal));
    deltaomega = omega(2) - omega(1);
    kimag(omega==omega1) = kimag1;
    
    % % Primer punto (excluye omega(1))
    for k = 2:g
        b(1) = b(1) + (nreal(k)-nreal(omega1)) / ((omega(k)^2 - omega(1)^2)*(omega(k)^2 - omega1^2));
    end
    kimag(1) = -2 * omega(1) / pi * deltaomega * b(1);

    % Último punto (excluye omega(g))
    for k = 1:g-1
        a(g) = a(g) + (nreal(k)-1) / (omega(k)^2 - omega(g)^2);
    end
    kimag(g) = -2 * omega(g) / pi * deltaomega * a(g);

    % Puntos intermedios
    for j = 2:g-1
        for k = 1:j-1
            a(j) = a(j) + (nreal(k)-1) / (omega(k)^2 - omega(j)^2);
        end
        for k = j+1:g
            b(j) = b(j) + (nreal(k)-1)/ (omega(k)^2 - omega(j)^2);
        end
        kimag(j) = -2 * omega(j) / pi * deltaomega * (a(j) + b(j));
    end
end