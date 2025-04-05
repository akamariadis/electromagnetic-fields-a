clear;
clc;
close all;

N = 50;
h_samples = linspace(0.01, 1, 200);
re_sigma = linspace(0, 2, N);
im_sigma = linspace(-2, 2, N);
[Re, Im] = meshgrid(re_sigma, im_sigma);
sigma_eta = Re + 1i*Im;

max_absorption = zeros(N);
best_h = zeros(N);

for i = 1:N
    for j = 1:N
        s = sigma_eta(i,j);
        absorptions = zeros(size(h_samples));
        for k = 1:length(h_samples)
            absorptions(k) = compute_absorption(s, h_samples(k));
        end
        [max_absorption(i,j), idx] = max(absorptions);
        best_h(i,j) = h_samples(idx);
    end
end

% (Γ)
figure;
imagesc(re_sigma, im_sigma, max_absorption);
axis xy; colorbar;
xlabel('Re[\sigma\eta_0]'); ylabel('Im[\sigma\eta_0]');
title('Μέγιστη Απορρόφηση');

figure;
imagesc(re_sigma, im_sigma, best_h);
axis xy; colorbar;
xlabel('Re[\sigma\eta_0]'); ylabel('Im[\sigma\eta_0]');
title('Βέλτιστο h/\lambda');

% (Δ)
best_idx = find(max_absorption == max(max_absorption(:)), 1);
[mid_idx, ~] = find(abs(max_absorption - 0.5) == min(abs(max_absorption(:)-0.5)), 1);
[bad_idx, ~] = find(max_absorption == min(max_absorption(:)), 1);

[idx1, idx2] = ind2sub([N N], [bad_idx mid_idx best_idx]);

labels = {'Κακό', 'Μέτριο', 'Βέλτιστο'};
for m = 1:3
    i = idx1(m); j = idx2(m);
    s = sigma_eta(i,j);
    h = best_h(i,j);
    plot_fields(s, h, labels{m});
end

% (Ε)
robustness_map(sigma_eta(idx1(3), idx2(3)), best_h(idx1(3), idx2(3)), re_sigma, im_sigma);

% (Ζ)
lossy_medium_analysis(sigma_eta(idx1(3), idx2(3)), 0.5);
