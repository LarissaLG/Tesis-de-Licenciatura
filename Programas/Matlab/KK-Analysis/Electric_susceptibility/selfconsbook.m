function [refin, imfin] = selfconsbook(omega, rechi, imchi, N, mu)
% selfconsbook - Estimación auto-consistente de la susceptibilidad compleja
% utilizando las relaciones de Kramers-Kronig.
%
% Entradas:
%   omega - vector de frecuencias o energías (eje horizontal)
%   rechi - vector de la parte real inicial (estimación o medida)
%   imchi - vector de la parte imaginaria inicial (estimación o medida)
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
    if size(rechi,1) > size(rechi,2)
        rechi = rechi';
    end
    if size(imchi,1) > size(imchi,2)
        imchi = imchi';
    end

    % Inicializa variables auxiliares (copias)
    comodo1 = rechi;   % Estimación de la parte real
    comodo2 = imchi;   % Estimación de la parte imaginaria

    % Bucle auto-consistente
    for j = 1:N
        % Calcula nueva parte real a partir de Im[χ] anterior
        comodo1 = kkrebook(omega, comodo2, 0);
        % Mezcla con la estimación inicial
        comodo1 = mu * rechi + (1 - mu) * comodo1;

        % Calcula nueva parte imaginaria a partir de Re[χ] actualizada
        comodo2 = kkimbook(omega, comodo1, 0);
        % Mezcla con la estimación inicial
        comodo2 = mu * imchi + (1 - mu) * comodo2;
    end

    % Devuelve las estimaciones auto-consistentes finales
    refin = comodo1;
    imfin = comodo2;
end
