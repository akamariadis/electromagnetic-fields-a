clear all;
close all;
clc;

% Constants
n_air = 1;
lambda = linspace(400, 700, 100) * 1e-9;
theta = linspace(0, 89, 90);
h = 0.8 * 530e-9;

% Materials (with vectorized refractive index functions)
materials = {'Ag', 'Au', 'a-Si', 'c-Si'};
n_material = containers.Map();

% Define material dispersion (ensure array compatibility)
n_material('Ag') = @(lam) 0.05 + 3.3i * (500e-9./lam).^2;
n_material('Au') = @(lam) 0.15 + 2.5i * (500e-9./lam).^2;
n_material('a-Si') = @(lam) 3.8 + 0.01i * ones(size(lam));
n_material('c-Si') = @(lam) 4.0 + 0.005i * ones(size(lam));

% Main loop
for mat_idx = 1:length(materials)
    material = materials{mat_idx};
    n = n_material(material)(lambda);
    
    % Preallocate reflection matrices
    R_TE = zeros(length(lambda), length(theta));
    R_TM = zeros(length(lambda), length(theta));
    
    % Calculate R(λ, θ) for TE and TM
    for lam_idx = 1:length(lambda)
        lam = lambda(lam_idx);
        n_lam = n(lam_idx);
        
        for theta_idx = 1:length(theta)
            th = theta(theta_idx);
            
            % TE polarization (s-polarized)
            [r_te, ~] = multilayer_reflection(n_air, n_lam, n_air, th, lam, h, 'TE');
            R_TE(lam_idx, theta_idx) = abs(r_te)^2;
            
            % TM polarization (p-polarized)
            [r_tm, ~] = multilayer_reflection(n_air, n_lam, n_air, th, lam, h, 'TM');
            R_TM(lam_idx, theta_idx) = abs(r_tm)^2;
        end
    end
    
    % Plot TE reflection
    figure;
    contourf(theta, lambda*1e9, R_TE, 20, 'LineColor', 'none');
    colorbar;
    title(['TE Polarization: |r|^2, ', material, ', h = 0.8\lambda_0']);
    xlabel('Angle of Incidence \theta (°)');
    ylabel('Wavelength \lambda (nm)');
    
    % Plot TM reflection
    figure;
    contourf(theta, lambda*1e9, R_TM, 20, 'LineColor', 'none');
    colorbar;
    title(['TM Polarization: |r|^2, ', material, ', h = 0.8\lambda_0']);
    xlabel('Angle of Incidence \theta (°)');
    ylabel('Wavelength \lambda (nm)');
end

% Helper function: Transfer matrix method for single layer
function [r, t] = multilayer_reflection(n1, n2, n3, theta_deg, lambda, h, pol)
    theta1 = theta_deg * pi/180;
    k0 = 2*pi / lambda;
    
    % Tangential wavevector (conserved)
    kx = n1 * k0 * sin(theta1);
    
    % Normal wavevectors in each layer
    kz1 = sqrt((n1*k0)^2 - kx^2);
    kz2 = sqrt((n2*k0)^2 - kx^2);
    kz3 = sqrt((n3*k0)^2 - kx^2);
    
    % Fresnel coefficients
    if strcmp(pol, 'TE')
        r12 = (kz1 - kz2) / (kz1 + kz2);
        r23 = (kz2 - kz3) / (kz2 + kz3);
    else
        r12 = ( (n2^2)*kz1 - (n1^2)*kz2 ) / ( (n2^2)*kz1 + (n1^2)*kz2 );
        r23 = ( (n3^2)*kz2 - (n2^2)*kz3 ) / ( (n3^2)*kz2 + (n2^2)*kz3 );
    end
    
    % Phase shift in the layer
    phi = kz2 * h;
    
    % Total reflection coefficient (single interface)
    r = (r12 + r23 * exp(2i * phi)) / (1 + r12 * r23 * exp(2i * phi));
    t = 0; % Transmission not needed here
end
