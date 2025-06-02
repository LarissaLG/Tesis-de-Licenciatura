function parte_re = kkrebook_prueba(omega_datos, parte_imag_datos)
    % omega_datos: vector de frecuencias (datos)
    % parte_imag_datos: parte imaginaria evaluada en esas frecuencias

    syms w w_prime real  % Usamos w en lugar de omega para evitar conflicto

    % Interpolamos los datos como una función simbólica (vía función anónima)
    F_interp = griddedInterpolant(omega_datos, parte_imag_datos, 'spline');
    imag_fun = matlabFunction(F_interp, 'Vars', w_prime);  % función anónima para la parte imaginaria

    % Definimos el integrando de Kramers-Kronig
    integrando = (w_prime * imag_fun(w_prime)) / (w_prime^2 - w^2);

    % Definimos la integral con valor principal
    KK_integral(w) = (2/pi) * int(integrando, w_prime, omega_datos(1), omega_datos(end), 'PrincipalValue', true);

    % Evaluamos en los puntos originales
    parte_re = double(KK_integral(omega_datos));  % omega_datos es numérico
end
