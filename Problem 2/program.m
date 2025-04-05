% Parameters
eta0 = 377;
lambda = 1;
h_vals = linspace(0.01, 3, 200);
re_vals = linspace(-10, 10, 100);
im_vals = linspace(-10, 10, 100);

% Create grids
[ReS, ImS] = meshgrid(re_vals, im_vals);
Pmax = zeros(size(ReS));
h_opt = zeros(size(ReS));

% Loop over each point in the complex conductivity plane
for i = 1:numel(ReS)
    sigma = (ReS(i) + 1j*ImS(i)) / eta0;  % σ = (σ*η0)/η0 = original σ
    Pabs = zeros(size(h_vals));

    for j = 1:length(h_vals)
        h = h_vals(j);
        k = 2*pi / lambda;
        Zs = 1 ./ sigma;
        Zin = -1j * eta0 * cot(k*h) + Zs;
        Gamma = (Zin - eta0) ./ (Zin + eta0);
        T = 2 * eta0 ./ (Zin + eta0);
        Pabs(j) = 1 - abs(Gamma).^2 - abs(T).^2;
    end

    % Find max absorption and corresponding h/lambda
    [Pmax(i), idx] = max(real(Pabs));  % Only real absorption
    h_opt(i) = h_vals(idx);
end

% Reshape results
Pmax_plot = reshape(Pmax, size(ReS));
hopt_plot = reshape(h_opt, size(ReS));

% Plot 1
figure;
imagesc(re_vals, im_vals, Pmax_plot);
xlabel('Re[\sigma]\cdot\eta_0');
ylabel('Im[\sigma]\cdot\eta_0');
title('Maximum Absorbed Power');
colorbar;
axis xy;
set(gca, 'FontSize', 12);

% Plot 2
figure;
imagesc(re_vals, im_vals, hopt_plot);
xlabel('Re[\sigma]\cdot\eta_0');
ylabel('Im[\sigma]\cdot\eta_0');
title('Optimal h/\lambda');
colorbar;
axis xy;
set(gca, 'FontSize', 12);

% Constants
eta0 = 377;
lambda = 1;
k = 2*pi / lambda;
z_vals = linspace(-1, 2, 1000); % Από z = -λ έως 2λ
E0 = 1;

sigmas_eta0 = [0.2+0.1j, 1+0.5j, 2+1j];
labels = {'Κακό', 'Μέτριο', 'Βέλτιστο'};
colors = {'r', 'g', 'b'};

figure; hold on;
for i = 1:3
    sigma = sigmas_eta0(i) / eta0;
    Zs = 1 / sigma;

    h = 0.25; 
    Zin = -1j * eta0 * cot(k*h) + Zs;
    Gamma = (Zin - eta0) / (Zin + eta0);
    T = 2 * eta0 / (Zin + eta0);

    E = zeros(size(z_vals));
    for j = 1:length(z_vals)
        z = z_vals(j);
        if z < 0
            E(j) = E0 * (exp(-1j*k*z) + Gamma * exp(1j*k*z));
        else
            E(j) = T * E0 * exp(-1j*k*z);
        end
    end

    plot(z_vals, abs(E), 'Color', colors{i}, 'DisplayName', labels{i});
end

xlabel('z/\lambda');
ylabel('|E(z)|');
title('Μέτρο Ηλεκτρικού Πεδίου για 3 Σενάρια');
legend show;
grid on;
set(gca, 'FontSize', 12);

% Constants
eta0 = 377;
lambda = 1;
k = 2*pi / lambda;
sigma_opt = (2 + 1j) / eta0;
h_opt = 0.25;

re_range = linspace(1, 3, 100);   % γύρω από Re[σ]*η0 = 2
im_range = linspace(0, 2, 100);   % γύρω από Im[σ]*η0 = 1
[ReS, ImS] = meshgrid(re_range, im_range);

AbsorbMap = zeros(size(ReS));
for i = 1:numel(ReS)
    sigma = (ReS(i) + 1j*ImS(i)) / eta0;
    Zs = 1 / sigma;
    Zin = -1j * eta0 * cot(k * h_opt) + Zs;
    Gamma = (Zin - eta0) / (Zin + eta0);
    T = 2 * eta0 / (Zin + eta0);
    AbsorbMap(i) = 1 - abs(Gamma)^2 - abs(T)^2;
end

% Create Plot 1

figure;
imagesc(re_range, im_range, real(AbsorbMap));
xlabel('Re[\sigma]\cdot\eta_0');
ylabel('Im[\sigma]\cdot\eta_0');
title('Απορροφητικότητα γύρω από το βέλτιστο \sigma (σταθερό h/\lambda)');
colorbar;
axis xy;
set(gca, 'FontSize', 12);

h_vals = linspace(0.01, 1, 300);
Absorb_h = zeros(size(h_vals));
for i = 1:length(h_vals)
    h = h_vals(i);
    Zs = 1 / sigma_opt;
    Zin = -1j * eta0 * cot(k * h) + Zs;
    Gamma = (Zin - eta0) / (Zin + eta0);
    T = 2 * eta0 / (Zin + eta0);
    Absorb_h(i) = 1 - abs(Gamma)^2 - abs(T)^2;
end

% Create Plot 2

figure;
plot(h_vals, real(Absorb_h), 'b', 'LineWidth', 2);
xlabel('h/\lambda');
ylabel('Απορροφητικότητα');
title('Απορροφητικότητα ως προς h/\lambda για βέλτιστο \sigma');
grid on;
set(gca, 'FontSize', 12);

% Constants
eta0 = 377;
lambda = 1;
sigma = (2 + 1j) / eta0;
Zs = 1 / sigma;
h = 0.25;  % σταθερό h/λ

re_eps = linspace(1, 5, 100);
im_eps = linspace(0, 2, 100);
[ReEps, ImEps] = meshgrid(re_eps, im_eps);
AbsorbMap = zeros(size(ReEps));

for i = 1:numel(ReEps)
    eps_r = ReEps(i) + 1j * ImEps(i);
    k = 2*pi / lambda * sqrt(eps_r);
    eta = eta0 / sqrt(eps_r);
    Zin = -1j * eta * cot(k * h) + Zs + Zs;
    Gamma = (Zin - eta0) / (Zin + eta0);
    T = 2 * eta0 / (Zin + eta0);
    AbsorbMap(i) = 1 - abs(Gamma)^2 - abs(T)^2;
end

% Create heatmap
figure;
imagesc(re_eps, im_eps, real(AbsorbMap));
xlabel('Re[\epsilon_r]');
ylabel('Im[\epsilon_r]');
title('Απορροφητικότητα με ενδιάμεσο υλικό απωλειών (h/\lambda = 0.25)');
colorbar;
axis xy;
set(gca, 'FontSize', 12);
