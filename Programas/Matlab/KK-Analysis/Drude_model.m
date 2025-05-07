function [epsilon, n_complex] = drude_model(omega, eps_inf, omega_p, gamma)
% DRUDE_MODEL - Calcula la epsilon y el índice de refracción usando el modelo de Drude
%
% Inputs:
%   omega    -  (rad/s)
%   eps_inf  - High-frequency dielectric constant (unitless)
%   omega_p  - Plasma frequency (rad/s)
%   gamma    - Damping rate (rad/s)
%
% Outputs:
%   epsilon     - Complex dielectric function epsilon(omega)
%   n_complex   - Complex refractive index: n + i*k

    % Drude model for complex dielectric function
    epsilon = eps_inf - (omega_p^2) ./ (omega.^2 + 1i * gamma .* omega);

    % Complex refractive index
    n_complex = sqrt(epsilon);
end

%%
% Constants
c = 3e8; % Speed of light (m/s)
hbar = 6.582119569 * 1e-16;

% Define wavelength range (in microns and meters)
lambda_um = linspace(0.05, 200, 500); % 0.2 µm to 2 µm
lambda_nm = lambda_um*1000; % 0.2 µm to 2 µm
lambda_m = lambda_um * 1e-6;

% Convert to angular frequency (rad/s)
omega = 2 * pi * c ./ lambda_m;


% Drude parameters (example values for gold)
eps_inf = 1;
omega_p = 13.142/hbar; % rad/s
gamma =  0.197/hbar;   % rad/s

[eps, n_complex] = drude_model(omega, eps_inf,omega_p, gamma);
n = real(n_complex);
k = imag(n_complex);

% Plot
% figure;
% plot(lambda_nm, n, 'b', 'LineWidth', 2); hold on;
% plot(lambda_nm, k, 'r', 'LineWidth', 2);
% xlabel('\lambda (\mum)');
% ylabel('n, k');
% legend('n(\lambda)', 'k(\lambda)');
% title('Drude Model: Complex Refractive Index');
% grid on;


%%
[refin, imfin]= selfconsbook_refractive_index(omega, n, k, 100,1);
%%
figure;
subplot(2,1,1)
plot(lambda_nm, n, 'b',lambda_nm, refin, 'r--')
xlabel('\lambda [nm]'); ylabel('n');
legend('n original','n estimado KK'); title('Parte real del índice');

subplot(2,1,2)
plot(lambda_nm, k, 'b',lambda_nm, imfin, 'r--')
xlabel('\lambda [nm]'); ylabel('k');
legend('k original','k estimado KK'); title('Parte imaginaria del índice');
