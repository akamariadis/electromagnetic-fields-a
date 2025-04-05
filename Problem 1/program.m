% Μήκη κύματος στο ορατό φάσμα (σε μικρόμετρα)
lambda_micron = linspace(0.4, 0.7, 100); 

% Υλικά που εξετάζουμε
metals = {'Ag', 'Au', 'Ti', 'Be', 'Al'};
dielectrics = {'aSi', 'cSi', 'Ge', 'GaAs', 'GaP'};
materials = [metals, dielectrics];

data = struct(); % Δομή αποθήκευσης δεδομένων

for i = 1:length(materials)
    % Φόρτωση δεδομένων από αρχείο (χρειάζεται να κατεβάσουμε τα αρχεία)
    filename = sprintf('%s_data.txt', materials{i});
    M = load(filename);
    
    lambda_data = M(:,1); % Μήκος κύματος (μm)
    n_data = M(:,2); % Πραγματικό μέρος
    k_data = M(:,3); % Φανταστικό μέρος
    
    % Παρεμβολή δεδομένων για το επιθυμητό εύρος
    n_interp = interp1(lambda_data, n_data, lambda_micron, 'linear', 'extrap');
    k_interp = interp1(lambda_data, k_data, lambda_micron, 'linear', 'extrap');
    
    % Υπολογισμός μιγαδικής επιτρεπτότητας
    epsilon = (n_interp - 1i * k_interp).^2;
    
    % Αποθήκευση
    data.(materials{i}) = struct('lambda', lambda_micron, 'n', n_interp, 'k', k_interp, 'epsilon', epsilon);
end

% Σχεδίαση για μέταλλα
figure;
hold on;
for i = 1:length(metals)
    plot(real(data.(metals{i}).epsilon), imag(data.(metals{i}).epsilon), 'LineWidth', 2);
end
xlabel('Re(\epsilon)');
ylabel('Im(\epsilon)');
title('Γεωμετρικοί τόποι επιτρεπτοτήτων - Μέταλλα');
legend(metals);
grid on;
hold off;

% Σχεδίαση για διηλεκτρικά
figure;
hold on;
for i = 1:length(dielectrics)
    plot(real(data.(dielectrics{i}).epsilon), imag(data.(dielectrics{i}).epsilon), 'LineWidth', 2);
end
xlabel('Re(\epsilon)');
ylabel('Im(\epsilon)');
title('Γεωμετρικοί τόποι επιτρεπτοτήτων - Διηλεκτρικά');
legend(dielectrics);
grid on;
hold off;
num_materials = length(metals);
optimal_ratios = zeros(num_materials);

for i = 1:num_materials
    for j = 1:num_materials
        eps1 = data.(metals{i}).epsilon;
        eps2 = data.(dielectrics{j}).epsilon;
        
        % Εύρεση f που κάνει την επιτρεπτότητα κοντά στο 1
        f_opt = (1 - eps2) ./ (eps1 - eps2);
        f_opt = max(min(f_opt, 1), 0); % Περιορισμός στο [0,1]
        
        optimal_ratios(i, j) = mean(f_opt); % Μέση τιμή
    end
end

% Εμφάνιση αποτελεσμάτων
figure;
imagesc(real(optimal_ratios)); % Παίρνουμε μόνο το πραγματικό μέρος
colorbar;
title('Βέλτιστες αναλογίες μετάλλου/διηλεκτρικού');
xlabel('Διηλεκτρικά'); ylabel('Μέταλλα');
% --- Επανάληψη του (Β) χωρίς απώλειες ---
optimal_ratios_no_loss = zeros(num_materials);

for i = 1:num_materials
    for j = 1:num_materials
        eps1 = real(data.(metals{i}).epsilon); % Πραγματικό μέρος (χωρίς απώλειες)
        eps2 = real(data.(dielectrics{j}).epsilon);
        
        % Εύρεση βέλτιστης αναλογίας f ώστε η επιτρεπτότητα να είναι κοντά στο 1
        f_opt = (1 - eps2) ./ (eps1 - eps2);
        f_opt = max(min(f_opt, 1), 0); % Περιορισμός στο [0,1]
        
        optimal_ratios_no_loss(i, j) = mean(f_opt); % Μέση τιμή
    end
end

% Σχεδίαση του πίνακα διαστρωματώσεων
figure;
imagesc(optimal_ratios_no_loss);
colorbar;
title('Βέλτιστες αναλογίες χωρίς απώλειες');
xlabel('Διηλεκτρικά'); ylabel('Μέταλλα');
Fs = 1e15; % Συχνότητα δειγματοληψίας
T = 1/Fs;
t = linspace(-1e-12, 1e-12, 1000); % Χρονικό εύρος

% --- Υπολογισμός κρουστικής απόκρισης για μέταλλα ---
figure; hold on;
for i = 1:length(metals)
    eps_w = data.(metals{i}).epsilon;
    h_eps_t = ifftshift(ifft(eps_w, length(t))); 

    h_eps_t_resized = interp1(linspace(-1e-12, 1e-12, length(h_eps_t)), h_eps_t, t, 'linear', 'extrap');

    plot(t, real(h_eps_t_resized), 'LineWidth', 2);
end
xlabel('Χρόνος (s)'); % Χωρίς LaTeX
ylabel('$h_{\varepsilon}(t)$', 'Interpreter', 'latex'); % LaTeX μόνο εδώ
title('Κρουστική απόκριση - Μέταλλα'); % Χωρίς LaTeX
legend(metals);
grid on;
hold off;

% --- Υπολογισμός κρουστικής απόκρισης για διηλεκτρικά ---
figure; hold on;
for i = 1:length(dielectrics)
    eps_w = data.(dielectrics{i}).epsilon;
    h_eps_t = ifftshift(ifft(eps_w, length(t)));

    h_eps_t_resized = interp1(linspace(-1e-12, 1e-12, length(h_eps_t)), h_eps_t, t, 'linear', 'extrap');

    plot(t, real(h_eps_t_resized), 'LineWidth', 2);
end
xlabel('Χρόνος (s)'); % Χωρίς LaTeX
ylabel('$h_{\varepsilon}(t)$', 'Interpreter', 'latex'); % LaTeX μόνο εδώ
title('Κρουστική απόκριση - Διηλεκτρικά'); % Χωρίς LaTeX
legend(dielectrics);
grid on;
hold off;
% Παράμετροι παλμού
omega0 = 2 * pi * 3e8 / (530e-9);
Omega = omega0 / 20;
t = linspace(-1e-12, 1e-12, 1000); % Χρονικό εύρος

% Gaussian παλμός
E_t = exp(-(Omega * t).^2) .* cos(omega0 * t);

% --- Απόκριση για μέταλλα ---
figure; hold on;
for i = 1:length(metals)
    eps_t = ifftshift(ifft(data.(metals{i}).epsilon, length(t))); % IFFT με σωστό μήκος
    D_t = eps_t .* E_t;

    plot(t, real(D_t), 'LineWidth', 2);
end
xlabel('Χρόνος (s)');
ylabel('D(t)/\epsilon_0');
title('Απόκριση Gaussian Παλμού - Μέταλλα');
legend(metals);
grid on;
hold off;