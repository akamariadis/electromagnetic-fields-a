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

clear; clc;

lambda0 = 530e-9;
h = 1.6 * lambda0;
theta1 = deg2rad(45);
n1 = 1.0;
n3 = 1.5;

materials = {'Ag', 0.18 + 2.95i, 'metal'; ...
             'Ti', 2.78 + 3.23i, 'metal'; ...
             'GaP', 3.3 + 0i, 'dielectric'; ...
             'c-Si', 3.88 + 0.02i, 'dielectric'};

z = linspace(-lambda0, h + lambda0, 1000);
k0 = 2*pi / lambda0;

E_fields = zeros(length(materials), length(z));

for m = 1:length(materials)
    n2 = materials{m,2};

    theta2 = asin(n1 * sin(theta1) / n2);
    theta3 = asin(n1 * sin(theta1) / n3);

    k1z = k0 * n1 * cos(theta1);
    k2z = k0 * n2 * cos(theta2);
    k3z = k0 * n3 * cos(theta3);

    r12 = (n2^2 * cos(theta1) - n1^2 * cos(theta2)) / (n2^2 * cos(theta1) + n1^2 * cos(theta2));
    r23 = (n3^2 * cos(theta2) - n2^2 * cos(theta3)) / (n3^2 * cos(theta2) + n2^2 * cos(theta3));

    delta = k2z * h;

    r_tot = (r12 + r23 * exp(2i * delta)) / (1 + r12 * r23 * exp(2i * delta));
    t01 = 1 + r_tot;

    for i = 1:length(z)
        zi = z(i);
        if zi < 0
            E_fields(m,i) = exp(1i * k1z * zi) + r_tot * exp(-1i * k1z * zi);
        elseif zi <= h
            E_fields(m,i) = t01 * (exp(1i * k2z * zi) + r23 * exp(1i * k2z * (2*h - zi)));
        else
            t12 = 1 + r23;
            E_fields(m,i) = t01 * t12 * exp(1i * k3z * (zi - h));
        end
    end
end

figure;
colors = {'b', 'r--', 'g', 'm--'};
hold on;
for m = 1:length(materials)
    plot(z*1e9, abs(E_fields(m,:)), colors{m}, 'LineWidth', 1.5);
end
xline(0, '--k');
xline(h*1e9, '--k');
xlabel('z (nm)');
ylabel('|E(z)|');
title('|E(z)| για TM Πόλωση – Μέταλλα & Διηλεκτρικά', 'FontWeight', 'bold');
legend('Ag (καλύτερο μέταλλο)', 'Ti (χειρότερο μέταλλο)', ...
       'GaP (καλύτερο διηλεκτρικό)', 'c-Si (χειρότερο διηλεκτρικό)');
grid on;

clear; clc;

lambda0 = 530e-9;
h = 1.6 * lambda0;
theta1 = deg2rad(45);
n1 = 1.0;
n3 = 1.5;

materials = {'Ag', 0.18 + 2.95i, 'metal'; ...
             'Ti', 2.78 + 3.23i, 'metal'; ...
             'GaP', 3.3 + 0i, 'dielectric'; ...
             'c-Si', 3.88 + 0.02i, 'dielectric'};

z = linspace(-lambda0, h + lambda0, 1000);
k0 = 2*pi / lambda0;

E_fields = zeros(length(materials), length(z));

for m = 1:length(materials)
    n2 = materials{m,2};

    theta2 = asin(n1 * sin(theta1) / n2);
    theta3 = asin(n1 * sin(theta1) / n3);

    k1z = k0 * n1 * cos(theta1);
    k2z = k0 * n2 * cos(theta2);
    k3z = k0 * n3 * cos(theta3);

    r12 = (n1 * cos(theta1) - n2 * cos(theta2)) / (n1 * cos(theta1) + n2 * cos(theta2));
    r23 = (n2 * cos(theta2) - n3 * cos(theta3)) / (n2 * cos(theta2) + n3 * cos(theta3));

    delta = k2z * h;

    r_tot = (r12 + r23 * exp(2i * delta)) / (1 + r12 * r23 * exp(2i * delta));
    t01 = 1 + r_tot;

    for i = 1:length(z)
        zi = z(i);
        if zi < 0
            E_fields(m,i) = exp(1i * k1z * zi) + r_tot * exp(-1i * k1z * zi);
        elseif zi <= h
            E_fields(m,i) = t01 * (exp(1i * k2z * zi) + r23 * exp(1i * k2z * (2*h - zi)));
        else
            t12 = 1 + r23;
            E_fields(m,i) = t01 * t12 * exp(1i * k3z * (zi - h));
        end
    end
end

figure;
colors = {'b', 'r--', 'g', 'm--'};
hold on;
for m = 1:length(materials)
    plot(z*1e9, abs(E_fields(m,:)), colors{m}, 'LineWidth', 1.5);
end
xline(0, '--k');
xline(h*1e9, '--k');
xlabel('z (nm)');
ylabel('|E(z)|');
title('|E(z)| για TE Πόλωση – Μέταλλα & Διηλεκτρικά', 'FontWeight', 'bold');
legend('Ag (καλύτερο μέταλλο)', 'Ti (χειρότερο μέταλλο)', ...
       'GaP (καλύτερο διηλεκτρικό)', 'c-Si (χειρότερο διηλεκτρικό)');
grid on;
