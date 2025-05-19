clear; clc;

lambda_range = linspace(400e-9, 700e-9, 300);
theta_range = linspace(0, pi/2, 300);
[Lambda, Theta] = meshgrid(lambda_range, theta_range);
h = 0.8 * 530e-9;

lambda_data = [400 450 500 550 600 650 700]*1e-9;
n_data = [0.13 0.15 0.17 0.18 0.21 0.24 0.27];
k_data = [3.98 3.55 3.15 2.95 2.70 2.50 2.30];

n2_complex = interp1(lambda_data, n_data, Lambda) + ...
             1i*interp1(lambda_data, k_data, Lambda);

n1 = 1.0;
n3 = 1.5;

theta2 = asin(n1 * sin(Theta) ./ n2_complex);
theta3 = asin(n1 * sin(Theta) / n3);

r12 = (n2_complex.*cos(Theta) - n1*cos(theta2)) ./ ...
      (n2_complex.*cos(Theta) + n1*cos(theta2));
r23 = (n3*cos(theta2) - n2_complex.*cos(theta3)) ./ ...
      (n3*cos(theta2) + n2_complex.*cos(theta3));

delta = 2 * pi .* n2_complex .* h .* cos(theta2) ./ Lambda;

r_total = (r12 + r23 .* exp(2i*delta)) ./ ...
          (1 + r12 .* r23 .* exp(2i*delta));
R = abs(r_total).^2;

figure;
contourf(rad2deg(Theta), Lambda*1e9, R, 50, 'LineColor', 'none');
xlabel('Γωνία πρόσπτωσης (°)', 'FontSize', 12);
ylabel('Μήκος κύματος (nm)', 'FontSize', 12);
title('|R| για Ag - TM πόλωση', 'FontWeight', 'bold');
colorbar;
colormap(jet);

clear; clc;

lambda0 = 530e-9;
h = 1.6 * lambda0;
theta1 = deg2rad(45);
n1 = 1.0;
n3 = 1.5;

n_Ag = 0.18 + 2.95i;
n_Ti = 2.78 + 3.23i;

materials = {'Ag', n_Ag; 'Ti', n_Ti};

z = linspace(-lambda0, 2*lambda0, 1000);
E_fields = zeros(length(materials), length(z));

for m = 1:2
    name = materials{m, 1};
    n2 = materials{m, 2};

    theta2 = asin(n1 * sin(theta1) / n2);
    theta3 = asin(n1 * sin(theta1) / n3);

    k0 = 2*pi / lambda0;
    k1z = k0 * n1 * cos(theta1);
    k2z = k0 * n2 * cos(theta2);
    k3z = k0 * n3 * cos(theta3);

    r12 = (n2^2 * cos(theta1) - n1^2 * cos(theta2)) / ...
          (n2^2 * cos(theta1) + n1^2 * cos(theta2));
    r23 = (n3^2 * cos(theta2) - n2^2 * cos(theta3)) / ...
          (n3^2 * cos(theta2) + n2^2 * cos(theta3));

    t12 = 2*n2^2*cos(theta1) / ...
          (n2^2*cos(theta1) + n1^2*cos(theta2));

    delta = k2z * h;
    M11 = cos(delta);
    M12 = 1i * sin(delta) / (n2 * cos(theta2));
    M21 = 1i * n2 * cos(theta2) * sin(delta);
    M22 = cos(delta);

    r = (r12 + r23*exp(2i*delta)) / ...
        (1 + r12*r23*exp(2i*delta));

    t = t12 / (cos(delta) + 1i * (r23 - r12)*sin(delta)/(2));

    A = 1;
    B = r;

    for i = 1:length(z)
        zi = z(i);
        if zi < 0
            E_fields(m, i) = exp(-1i*k1z*zi) + r*exp(1i*k1z*zi);
        elseif zi <= h
            E_fields(m, i) = A*exp(-1i*k2z*zi) + B*exp(1i*k2z*zi);
        else
            E_fields(m, i) = t * exp(-1i*k3z*(zi - h));
        end
    end
end

figure;
plot(z*1e9, abs(E_fields(1,:)), 'b', 'LineWidth', 1.5); hold on;
plot(z*1e9, abs(E_fields(2,:)), 'r--', 'LineWidth', 1.5);
xlabel('z (nm)');
ylabel('|E(z)|');
legend('Ag (καλύτερο)', 'Ti (χειρότερο)');
title('|E(z)| για TM πόλωση - Μέταλλα', 'FontWeight', 'bold');
grid on;

clear; clc;

lambda = linspace(400e-9, 700e-9, 200); 
theta_deg = linspace(0, 80, 200);
[LL, TT] = meshgrid(lambda, theta_deg);

n1 = 1.0;
n3 = 1.5;

lambda0 = 530e-9;
h = 0.8 * lambda0;

n_Ag_data = @(l) 0.18 + 2.95i;
R = zeros(size(LL));

for i = 1:numel(LL)
    lam = LL(i);
    theta1 = deg2rad(TT(i));
    n2 = n_Ag_data(lam);

    theta2 = asin(n1*sin(theta1)/n2);
    theta3 = asin(n1*sin(theta1)/n3);

    r12 = (n1*cos(theta1) - n2*cos(theta2)) / ...
          (n1*cos(theta1) + n2*cos(theta2));
    r23 = (n2*cos(theta2) - n3*cos(theta3)) / ...
          (n2*cos(theta2) + n3*cos(theta3));

    delta = 2*pi*n2*cos(theta2)*h / lam;

    r_tot = (r12 + r23*exp(2i*delta)) / ...
            (1 + r12*r23*exp(2i*delta));

    R(i) = abs(r_tot)^2;
end

figure;
contourf(lambda*1e9, theta_deg, R, 50, 'LineColor', 'none');
xlabel('Μήκος κύματος λ (nm)');
ylabel('Γωνία πρόσπτωσης θ (°)');
title('|R(λ,θ)| για TE Πόλωση - Ag (h = 0.8λ₀)', 'FontWeight', 'bold');
colorbar;
colormap hot;

clear; clc;

lambda = linspace(400e-9, 700e-9, 200);
theta_deg = linspace(0, 80, 200);
[LL, TT] = meshgrid(lambda, theta_deg);

n1 = 1.0;
n3 = 1.5;

lambda0 = 530e-9;
h = 0.8 * lambda0;

n_GaP = 3.33 + 0i;   % μπορεί να γίνει συνάρτηση με dispersion

R = zeros(size(LL));

for i = 1:numel(LL)
    lam = LL(i);
    theta1 = deg2rad(TT(i));
    n2 = n_GaP;

    theta2 = asin(n1*sin(theta1)/n2);
    theta3 = asin(n1*sin(theta1)/n3);

    r12 = (n1*cos(theta1) - n2*cos(theta2)) / ...
          (n1*cos(theta1) + n2*cos(theta2));
    r23 = (n2*cos(theta2) - n3*cos(theta3)) / ...
          (n2*cos(theta2) + n3*cos(theta3));

    delta = 2*pi*n2*cos(theta2)*h / lam;

    r_tot = (r12 + r23*exp(2i*delta)) / ...
            (1 + r12*r23*exp(2i*delta));

    R(i) = abs(r_tot)^2;
end

figure;
contourf(lambda*1e9, theta_deg, R, 50, 'LineColor', 'none');
xlabel('Μήκος κύματος λ (nm)');
ylabel('Γωνία πρόσπτωσης θ (°)');
title('|R(λ,θ)| για TE Πόλωση - GaP (h = 0.8λ₀)', 'FontWeight', 'bold');
colorbar;
colormap turbo;

clear; clc;

lambda = linspace(400e-9, 700e-9, 200);
theta_deg = linspace(0, 80, 200);
[LL, TT] = meshgrid(lambda, theta_deg);

n1 = 1.0;
n3 = 1.5;

lambda0 = 530e-9;
h = 0.8 * lambda0;

n_Au_data = @(l) 0.5 + 2.9i;

R = zeros(size(LL));

for i = 1:numel(LL)
    lam = LL(i);
    theta1 = deg2rad(TT(i));
    n2 = n_Au_data(lam);

    theta2 = asin(n1*sin(theta1)/n2);
    theta3 = asin(n1*sin(theta1)/n3);

    r12 = (n1*cos(theta1) - n2*cos(theta2)) / ...
          (n1*cos(theta1) + n2*cos(theta2));
    r23 = (n2*cos(theta2) - n3*cos(theta3)) / ...
          (n2*cos(theta2) + n3*cos(theta3));

    delta = 2*pi*n2*cos(theta2)*h / lam;

    r_tot = (r12 + r23*exp(2i*delta)) / ...
            (1 + r12*r23*exp(2i*delta));

    R(i) = abs(r_tot)^2;
end

figure;
contourf(lambda*1e9, theta_deg, R, 50, 'LineColor', 'none');
xlabel('Μήκος κύματος λ (nm)');
ylabel('Γωνία πρόσπτωσης θ (°)');
title('|R(λ,θ)| για TE Πόλωση - Χρυσός (Au, h = 0.8λ₀)', 'FontWeight', 'bold');
colorbar;
colormap parula;

clear; clc;

lambda = linspace(400e-9, 700e-9, 200);
theta_deg = linspace(0, 80, 200);
[LL, TT] = meshgrid(lambda, theta_deg);

n1 = 1.0;
n3 = 1.5;

lambda0 = 530e-9;
h = 0.8 * lambda0;

n_cSi_data = @(l) 3.88 + 0.02i;

R = zeros(size(LL));

for i = 1:numel(LL)
    lam = LL(i);
    theta1 = deg2rad(TT(i));
    n2 = n_cSi_data(lam);

    theta2 = asin(n1*sin(theta1)/n2);
    theta3 = asin(n1*sin(theta1)/n3);

    r12 = (n1*cos(theta1) - n2*cos(theta2)) / ...
          (n1*cos(theta1) + n2*cos(theta2));
    r23 = (n2*cos(theta2) - n3*cos(theta3)) / ...
          (n2*cos(theta2) + n3*cos(theta3));

    delta = 2*pi*n2*cos(theta2)*h / lam;

    r_tot = (r12 + r23*exp(2i*delta)) / ...
            (1 + r12*r23*exp(2i*delta));

    R(i) = abs(r_tot)^2;
end

figure;
contourf(lambda*1e9, theta_deg, R, 50, 'LineColor', 'none');
xlabel('Μήκος κύματος λ (nm)');
ylabel('Γωνία πρόσπτωσης θ (°)');
title('|R(λ,θ)| για TE Πόλωση - c-Si (h = 0.8λ₀)', 'FontWeight', 'bold');
colorbar;
colormap autumn;

clear; clc;

lambda0 = 530e-9;
h = 1.6 * lambda0;
z = linspace(-lambda0, h + lambda0, 1000);
k0 = 2*pi/lambda0;

n1 = 1.0;
n3 = 1.5;

n_GaP = 3.3 + 0i;
n_cSi = 3.88 + 0.02i;

compute_Ez = @(n2) compute_field(z, lambda0, h, n1, n2, n3, k0);

E_GaP = compute_Ez(n_GaP);
E_cSi = compute_Ez(n_cSi);

figure;
plot(z*1e9, abs(E_GaP), 'b', 'LineWidth', 1.5); hold on;
plot(z*1e9, abs(E_cSi), 'r', 'LineWidth', 1.5);
xline(0, '--k'); xline(h*1e9, '--k');
xlabel('z (nm)');
ylabel('|E(z)|');
title('|E(z)| για GaP και c-Si (TE Πόλωση, λ = 530nm, h = 1.6λ₀)');
legend('GaP', 'c-Si');
grid on;

clear; clc;

lambda = linspace(400e-9, 700e-9, 150);
theta_deg = linspace(0, 80, 150);
[LL, TT] = meshgrid(lambda, theta_deg);

n1 = 1.0;
n3 = 1.5;
lambda0 = 530e-9;
h = 0.8 * lambda0;

materials = struct(...
    'name', {'Ag', 'Au', 'GaP', 'cSi'}, ...
    'n', {0.13 + 3.98i, 0.5 + 2.9i, 3.3 + 0i, 3.88 + 0.02i});

figure;
for k = 1:length(materials)
    R = zeros(size(LL));
    n2 = materials(k).n;

    for i = 1:numel(LL)
        lam = LL(i);
        theta1 = deg2rad(TT(i));

        % Snell
        theta2 = asin(n1*sin(theta1)/n2);
        theta3 = asin(n1*sin(theta1)/n3);

        % Fresnel TM
        r12 = (n2*cos(theta1) - n1*cos(theta2)) / (n2*cos(theta1) + n1*cos(theta2));
        r23 = (n3*cos(theta2) - n2*cos(theta3)) / (n3*cos(theta2) + n2*cos(theta3));

        delta = 2*pi*n2*cos(theta2)*h / lam;

        r_tot = (r12 + r23*exp(2i*delta)) / (1 + r12*r23*exp(2i*delta));

        R(i) = abs(r_tot)^2;
    end

    subplot(2,2,k)
    contourf(lambda*1e9, theta_deg, R, 50, 'LineColor', 'none');
    title(['|R(λ,θ)| TM - ' materials(k).name ' (h=0.8λ_0)']);
    xlabel('Μήκος κύματος λ (nm)');
    ylabel('Γωνία πρόσπτωσης θ (°)');
    colorbar;
    colormap turbo;
end

sgtitle('Contour plots συντελεστή ανάκλασης για TM πόλωση');

clear; clc;

lambda = linspace(400e-9, 700e-9, 150);
theta_deg = linspace(0, 80, 150);
[LL, TT] = meshgrid(lambda, theta_deg);

n1 = 1.0;
n3 = 1.5;
lambda0 = 530e-9;
h = 0.8 * lambda0;

materials = struct(...
    'name', {'Ag', 'Au', 'GaP', 'cSi'}, ...
    'n', {0.13 + 3.98i, 0.5 + 2.9i, 3.3 + 0i, 3.88 + 0.02i});

figure;
for k = 1:length(materials)
    R = zeros(size(LL));
    n2 = materials(k).n;

    for i = 1:numel(LL)
        lam = LL(i);
        theta1 = deg2rad(TT(i));

        theta2 = asin(n1*sin(theta1)/n2);
        theta3 = asin(n1*sin(theta1)/n3);

        r12 = (n1*cos(theta1) - n2*cos(theta2)) / (n1*cos(theta1) + n2*cos(theta2));
        r23 = (n2*cos(theta2) - n3*cos(theta3)) / (n2*cos(theta2) + n3*cos(theta3));

        delta = 2*pi*n2*cos(theta2)*h / lam;

        r_tot = (r12 + r23*exp(2i*delta)) / (1 + r12*r23*exp(2i*delta));

        R(i) = abs(r_tot)^2;
    end

    subplot(2,2,k)
    contourf(lambda*1e9, theta_deg, R, 50, 'LineColor', 'none');
    title(['|R(λ,θ)| TE - ' materials(k).name ' (h=0.8λ_0)']);
    xlabel('Μήκος κύματος λ (nm)');
    ylabel('Γωνία πρόσπτωσης θ (°)');
    colorbar;
    colormap turbo;
end

sgtitle('Contour plots συντελεστή ανάκλασης για TE πόλωση');

clear; clc;

lambda0 = 530e-9;
h = 1.6 * lambda0;
z = linspace(-lambda0, h + lambda0, 1000);

n1 = 1.0;
n3 = 1.5;

n_Ag = 0.13 + 3.98i;
n_Ti = 2.9 + 3.4i;

E_TM = @(z, n2) compute_E_TM(z, lambda0, h, n1, n2, n3);

E_Ag = E_TM(z, n_Ag);
E_Ti = E_TM(z, n_Ti);

figure;
plot(z*1e9, abs(E_Ag), 'b', 'LineWidth',1.5); hold on;
plot(z*1e9, abs(E_Ti), 'r', 'LineWidth',1.5);
xline(0, '--k');
xline(h*1e9, '--k');
xlabel('z (nm)');
ylabel('|E(z)|');
title('Κατανομή |E(z)| TM Πόλωση - Άργυρος (μπλε) & Τιτάνιο (κόκκινο)');
legend('Ag (καλύτερο)', 'Ti (χειρότερο)');
grid on;

n_GaP = 3.3 + 0i;
n_cSi = 3.88 + 0.02i;

E_GaP = E_TM(z, n_GaP);
E_cSi = E_TM(z, n_cSi);

figure;
plot(z*1e9, abs(E_GaP), 'b', 'LineWidth',1.5); hold on;
plot(z*1e9, abs(E_cSi), 'r', 'LineWidth',1.5);
xline(0, '--k');
xline(h*1e9, '--k');
xlabel('z (nm)');
ylabel('|E(z)|');
title('Κατανομή |E(z)| TM Πόλωση - GaP (μπλε) & c-Si (κόκκινο)');
legend('GaP (καλύτερο)', 'c-Si (χειρότερο)');
grid on;


% --- Συνάρτηση υπολογισμού πεδίου TM ---
function E = compute_E_TM(z, lambda, h, n1, n2, n3)
    k0 = 2*pi/lambda;
    theta1 = 0;

    theta2 = asin(n1*sin(theta1)/n2);
    theta3 = asin(n1*sin(theta1)/n3);

    r12 = (n2*cos(theta1) - n1*cos(theta2)) / (n2*cos(theta1) + n1*cos(theta2));
    r23 = (n3*cos(theta2) - n2*cos(theta3)) / (n3*cos(theta2) + n2*cos(theta3));
    delta = k0*n2*cos(theta2)*h;

    r_tot = (r12 + r23*exp(2i*delta)) / (1 + r12*r23*exp(2i*delta));
    t01 = 1 + r_tot;

    E = zeros(size(z));
    % z < 0 (αέρας)
    ind_air = z < 0;
    E(ind_air) = exp(1i*k0*n1*cos(theta1)*z(ind_air)) + r_tot*exp(-1i*k0*n1*cos(theta1)*z(ind_air));
    % 0 <= z <= h (στρώμα)
    ind_layer = (z >= 0) & (z <= h);
    beta2 = k0*n2*cos(theta2);
    E(ind_layer) = t01 * (exp(1i*beta2*z(ind_layer)) + r23*exp(1i*beta2*(2*h - z(ind_layer))));
    % z > h (υπόστρωμα)
    ind_sub = z > h;
    t12 = (1 + r23);
    beta3 = k0*n3*cos(theta3);
    E(ind_sub) = t01 * t12 * exp(1i*beta3*(z(ind_sub) - h));
end
