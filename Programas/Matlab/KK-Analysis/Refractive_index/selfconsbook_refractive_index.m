function [refin, imfin] = selfconsbook_refractive_index(omega, reN, imN, N,mu)
% selfconsbook - Estimación auto-consistente del índice de refracción
% complejo utilizando las relaciones de Kramers-Kronig.
%
% Entradas:
%   omega - vector de frecuencias o energías (eje horizontal)
%   reN - vector de la parte real inicial (estimación o medida)
%   imN - vector de la parte imaginaria inicial (estimación o medida)
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
    if size(reN,1) > size(reN,2)
        reN = reN';
    end
    if size(imN,1) > size(imN,2)
        imN = imN';
    end

    % Inicializa variables auxiliares (copias)
    comodo1 = reN;   % Estimación de la parte real
    comodo2 = imN;   % Estimación de la parte imaginaria
    max_iter = 100;  % Número máximo de iteraciones
    tol = 1e-6;  % tolerancia de convergencia

    % for j = 1:max_iter
    %     old1 = comodo1;
    %     old2 = comodo2;
    % 
    %     % Actualiza parte real
    %     comodo1 = kkrebook_refractive_index(omega, comodo2);
    %     comodo1 = mu * reN + (1 - mu) * comodo1;
    % 
    %     % Actualiza parte imaginaria
    %     comodo2 = kkimbook_refractive_index(omega, comodo1);
    %     comodo2 = mu * imN + (1 - mu) * comodo2;
    % 
    %     % Verifica convergencia relativa
    %     err1 = norm(comodo1 - old1) / norm(old1);
    %     err2 = norm(comodo2 - old2) / norm(old2);
    % 
    %     if err1 < tol && err2 < tol
    %         fprintf('Convergencia en %d iteraciones.\n', j);
    %         break
    %     end
    % end



    % Bucle auto-consistente
    for j = 1:N
        % Calcula nueva parte real a partir de Re[N] anterior
        comodo1 = kkrebook_refractive_index(omega, comodo2);
        % Mezcla con la estimación inicial
        comodo1 = mu * reN + (1 - mu) * comodo1;

        % Calcula nueva parte imaginaria a partir de Im[N] actualizada
        comodo2 = kkimbook_refractive_index(omega, comodo1);
        % Mezcla con la estimación inicial
        comodo2 = mu * imN + (1 - mu) * comodo2;
    end

    % Devuelve las estimaciones auto-consistentes finales
    refin = comodo1;
    imfin = comodo2;
end
