function [refin, imfin] = selfconsbook(omega, nreal, kimag, N,mu)
% selfconsbook - Estimación auto-consistente del índice de refracción
% complejo utilizando las relaciones de Kramers-Kronig.
%
% Entradas:
%   omega - vector de frecuencias o energías (eje horizontal)
%   nreal - vector de la parte real inicial (estimación o medida)
%   kimag - vector de la parte imaginaria inicial (estimación o medida)
%   N     - número de iteraciones del proceso auto-consistente
%   mu    - factor de mezcla entre el valor inicial y el nuevo (0 < mu < 1)
%
% Salidas:
%   refin - parte real auto-consistente de la susceptibilidad
%   imfin - parte imaginaria auto-consistente de la susceptibilidad

    % Asegura que todos los vectores sean filas
    if size(omega,1) > size(omega,2)
        omega = omega';
    end
    if size(nreal,1) > size(nreal,2)
        nreal = nreal';
    end
    if size(kimag,1) > size(kimag,2)
        kimag = kimag';
    end

    % Inicializa variables auxiliares (copias)
    comodo1 = nreal;   % Estimación de la parte real
    comodo2 = kimag;   % Estimación de la parte imaginaria

    % Bucle auto-consistente
    for j = 1:N
        % Calcula nueva parte real a partir de Re[N] anterior
        comodo1 = kkrebook(omega, comodo2);
        % Mezcla con la estimación inicial
        comodo1 = mu * nreal + (1 - mu) * comodo1;

        % Calcula nueva parte imaginaria a partir de Im[N] actualizada
        comodo2 = kkimbook(omega, comodo1);
        % Mezcla con la estimación inicial
        comodo2 = mu * kimag + (1 - mu) * comodo2;
    end

    % Devuelve las estimaciones auto-consistentes finales
    refin = comodo1;
    imfin = comodo2;
end
