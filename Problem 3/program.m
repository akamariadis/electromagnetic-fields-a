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
