clear; clc;

lambda_range = linspace(400e-9, 700e-9, 300); % [m]
theta_range = linspace(0, pi/2, 300);         % [rad]
[Lambda, Theta] = meshgrid(lambda_range, theta_range);
h = 0.8 * 530e-9; % πάχος στρώματος

lambda_data = [400 450 500 550 600 650 700]*1e-9;
n_data = [0.13 0.15 0.17 0.18 0.21 0.24 0.27];    % Re(n)
k_data = [3.98 3.55 3.15 2.95 2.70 2.50 2.30];    % Im(n)

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
