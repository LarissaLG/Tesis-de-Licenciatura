function kimag=sskkimbook(omega,nreal,omega1,kimag1)
% Calcula la parte imaginaria del índice de refracción/función dieléctrica
% a partir de la parte real n(omega), usando la relación sustractiva de KK.

% Entradas:
% omega - vector de energía en eVs (debe estar equiespaciado)
% nreal - vector de la parte real del índice de refracción/funcion
% dieléctrica
% omega1 - energía de kimag conocido
% kimag1 - valor de kimag conocido

% Salida:
% kimag - vector estimado de la parte imaginaria del índice de
% refracción/función dieléctrica

    if size(omega,1)>size(omega,2)
        omega = omega';
    end
    
    if size(nreal,1)>size(nreal,2)
        nreal = nreal';
    end
    
    g=size(omega,2);
 
    x = 0;

    % Encontrar el índice más cercano a omega1
    [~, x] = min(abs(omega - omega1));

    kimag = kkimbook(omega,nreal);
    %Aplicación de las relaciones de KK

    kimag = kimag + (omega / omega1) .* (kimag1 - kimag(x));

    %The subtracted relation upgrades the estimate obtained
    %with K-K relations.