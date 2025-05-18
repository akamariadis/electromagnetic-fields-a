clear all;
close all;
clc;

n_air = 1;
lambda = linspace(400, 700, 100) * 1e-9;
theta = linspace(0, 89, 90);
h = 0.8 * 530e-9;
materials = {'Ag', 'Au', 'a-Si', 'c-Si'};
n_material = containers.Map();
n_material('Ag') = @(lam) 0.05 + 3.3i * (500e-9./lam).^2;
n_material('Au') = @(lam) 0.15 + 2.5i * (500e-9./lam).^2;
n_material('a-Si') = @(lam) 3.8 + 0.01i * ones(size(lam));
n_material('c-Si') = @(lam) 4.0 + 0.005i * ones(size(lam));
for mat_idx = 1:length(materials)
    material = materials{mat_idx};
    n = n_material(material)(lambda);
    R_TE = zeros(length(lambda), length(theta));
    R_TM = zeros(length(lambda), length(theta));
    for lam_idx = 1:length(lambda)
        lam = lambda(lam_idx);
        n_lam = n(lam_idx);
        for theta_idx = 1:length(theta)
            th = theta(theta_idx);
            [r_te, ~] = multilayer_reflection(n_air, n_lam, n_air, th, lam, h, 'TE');
            R_TE(lam_idx, theta_idx) = abs(r_te)^2;
            [r_tm, ~] = multilayer_reflection(n_air, n_lam, n_air, th, lam, h, 'TM');
            R_TM(lam_idx, theta_idx) = abs(r_tm)^2;
        end
    end
    figure;
    contourf(theta, lambda*1e9, R_TE, 20, 'LineColor', 'none');
    colorbar;
    title(['TE Polarization: |r|^2, ', material, ', h = 0.8\lambda_0']);
    xlabel('Angle of Incidence \theta (°)');
    ylabel('Wavelength \lambda (nm)');
    figure;
    contourf(theta, lambda*1e9, R_TM, 20, 'LineColor', 'none');
    colorbar;
    title(['TM Polarization: |r|^2, ', material, ', h = 0.8\lambda_0']);
    xlabel('Angle of Incidence \theta (°)');
    ylabel('Wavelength \lambda (nm)');
end
function [r, t] = multilayer_reflection(n1, n2, n3, theta_deg, lambda, h, pol)
    theta1 = theta_deg * pi/180;
    k0 = 2*pi / lambda;
    kx = n1 * k0 * sin(theta1);
    kz1 = sqrt((n1*k0)^2 - kx^2);
    kz2 = sqrt((n2*k0)^2 - kx^2);
    kz3 = sqrt((n3*k0)^2 - kx^2);
    if strcmp(pol, 'TE')
        r12 = (kz1 - kz2) / (kz1 + kz2);
        r23 = (kz2 - kz3) / (kz2 + kz3);
    else
        r12 = ( (n2^2)*kz1 - (n1^2)*kz2 ) / ( (n2^2)*kz1 + (n1^2)*kz2 );
        r23 = ( (n3^2)*kz2 - (n2^2)*kz3 ) / ( (n3^2)*kz2 + (n2^2)*kz3 );
    end
    phi = kz2 * h;
    r = (r12 + r23 * exp(2i * phi)) / (1 + r12 * r23 * exp(2i * phi));
    t = 0; % Transmission not needed here
end
function plot_electric_field_profiles()
    lambda0 = 530e-9;
    h = 1.6 * lambda0;
    n_air = 1;
    materials_metal = {'Ag', 'Ti'};
    materials_diel = {'a-Si', 'c-Si'};
    polarizations = {'TE', 'TM'};
    lambda = lambda0;
    theta = 0;
    z = linspace(0, h, 1000);
    for pol_idx = 1:length(polarizations)
        pol = polarizations{pol_idx};
        figure;
        hold on;
        for mat_idx = 1:length(materials_metal)
            material = materials_metal{mat_idx};
            n = getRefractiveIndex(material, lambda);
            [E, ~] = calculate_field_profile(n_air, n, n_air, theta, lambda, h, z, pol);
            plot(z/lambda0, abs(E)/max(abs(E)), 'LineWidth', 2, 'DisplayName', material);
        end
        xlabel('Position z / \lambda_0');
        ylabel('|E(z)| / |E_{inc}|');
        title(['Metals: Electric Field Profile (', pol, ' Polarization)']);
        legend('Location', 'best');
        grid on
        figure;
        hold on;
        for mat_idx = 1:length(materials_diel)
            material = materials_diel{mat_idx};
            n = getRefractiveIndex(material, lambda);
            [E, ~] = calculate_field_profile(n_air, n, n_air, theta, lambda, h, z, pol);
            plot(z/lambda0, abs(E)/max(abs(E)), 'LineWidth', 2, 'DisplayName', material);
        end
        xlabel('Position z / \lambda_0');
        ylabel('|E(z)| / |E_{inc}|');
        title(['Dielectrics: Electric Field Profile (', pol, ' Polarization)']);
        legend('Location', 'best');
        grid on;
    end
end
function n = getRefractiveIndex(material, lambda)
    switch material
        case 'Ag'
            n = 0.05 + 3.3i * (530e-9/lambda);
        case 'Ti'
            n = 2.5 + 3.5i;
        case 'a-Si'
            n = 3.8 + 0.01i;
        case 'c-Si'
            n = 4.0 + 0.005i;
        otherwise
            error('Material not supported.');
    end
end
function [E, H] = calculate_field_profile(n1, n2, n3, theta_deg, lambda, h, z, pol)
    theta1 = theta_deg * pi/180;
    k0 = 2*pi / lambda;
    kx = n1 * k0 * sin(theta1);
    kz1 = sqrt((n1*k0)^2 - kx^2);
    kz2 = sqrt((n2*k0)^2 - kx^2);
    kz3 = sqrt((n3*k0)^2 - kx^2);
    if strcmp(pol, 'TE')
        r12 = (kz1 - kz2) / (kz1 + kz2);
        r23 = (kz2 - kz3) / (kz2 + kz3);
        t12 = 2 * kz1 / (kz1 + kz2);
        t23 = 2 * kz2 / (kz2 + kz3);
    else
        r12 = ( (n2^2)*kz1 - (n1^2)*kz2 ) / ( (n2^2)*kz1 + (n1^2)*kz2 );
        r23 = ( (n3^2)*kz2 - (n2^2)*kz3 ) / ( (n3^2)*kz2 + (n2^2)*kz3 );
        t12 = 2 * n1 * n2 * kz1 / ( (n2^2)*kz1 + (n1^2)*kz2 );
        t23 = 2 * n2 * n3 * kz2 / ( (n3^2)*kz2 + (n2^2)*kz3 );
    end
    phi = kz2 * h;
    r = (r12 + r23 * exp(2i * phi)) / (1 + r12 * r23 * exp(2i * phi));
    t = (t12 * t23 * exp(1i * phi)) / (1 + r12 * r23 * exp(2i * phi));
    E = zeros(size(z));
    air_region = (z < 0);
    E(air_region) = exp(1i * kz1 * z(air_region)) + r * exp(-1i * kz1 * z(air_region));
    layer_region = (z >= 0) & (z <= h);
    E(layer_region) = t12 * exp(1i * kz2 * z(layer_region)) + r23 * t12 * exp(1i * kz2 * (2*h - z(layer_region)));
    sub_region = (z > h);
    E(sub_region) = t * exp(1i * kz3 * (z(sub_region) - h));
    H = zeros(size(z));
end
