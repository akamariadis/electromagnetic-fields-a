function main()
    n_air = 1;
    lambda = linspace(400, 700, 100) * 1e-9;
    theta = linspace(0, 89, 90);
    h = 0.8 * 530e-9;
    
    materials = {'Ag', 'Au', 'a-Si', 'c-Si'};
    n_Ag = getRefractiveIndex('Ag', lambda);
    n_Au = getRefractiveIndex('Au', lambda);
    n_aSi = getRefractiveIndex('a-Si', lambda);
    n_cSi = getRefractiveIndex('c-Si', lambda);
    n_data = containers.Map();
    n_data('Ag') = n_Ag;
    n_data('Au') = n_Au;
    n_data('a-Si') = n_aSi;
    n_data('c-Si') = n_cSi;
    
    for mat_idx = 1:length(materials)
        material = materials{mat_idx};
        n = n_data(material);
        R_TE = zeros(length(lambda), length(theta));
        R_TM = zeros(length(lambda), length(theta));        
        for lam_idx = 1:length(lambda)
            lam = lambda(lam_idx);
            n_lam = n(lam_idx);
            for theta_idx = 1:length(theta)
                th = theta(theta_idx);
                r_te = calculateReflection(n_air, n_lam, n_air, th, lam, h, 'TE');
                R_TE(lam_idx, theta_idx) = abs(r_te)^2;
                r_tm = calculateReflection(n_air, n_lam, n_air, th, lam, h, 'TM');
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
end

function n = getRefractiveIndex(material, lambda)
    switch material
        case 'Ag'
            n = 0.05 + 3.3i * (500e-9./lambda).^2;
        case 'Au'
            n = 0.15 + 2.5i * (500e-9./lambda).^2;
        case 'a-Si'
            n = 3.8 + 0.01i * ones(size(lambda));
        case 'c-Si'
            n = 4.0 + 0.005i * ones(size(lambda));
        otherwise
            error('Material not supported.');
    end
end

function r = calculateReflection(n1, n2, n3, theta_deg, lambda, h, pol)
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
end
