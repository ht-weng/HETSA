function [phi, error, pacf, significance] = estimateARbyBurgs(x, order)
% [phi_i, s2, C]=burg(x,p), estimation of AR(order) parameters 
% using Burg's algorithm, Brockwell p 146
% x time series, order AR order
% phi_i estimated parameter vector, 
% error estimated variance
% Gamma estimetad covariance matrix of the estimated pacf:s
    if(size(x, 1) > size(x, 2))
        disp('Make sure that timeseries are stored as rows');
        phi = [];
        error = [];
        pacf = [];
        significance = [];         
        return;
    end
    sampleSize = length(x);
    % d = x' * x - (x(1)^2 - x(n)^2);
    u = x(sampleSize:-1:1);
    v = u;
    d = 0;
    for k = 2:sampleSize
        d = d + (u(k-1).^2) + (v(k).^2);
    end
 
    for k = 1:order        
            % phi_kk(k) = v(2:n-k+1)' * u(1:n-k) / d;
            phi_kk(k) = ((v(k+1:sampleSize) * u(k:sampleSize-1)') * 2) / d;
    
            u1 = u;
            % forward one-step prediction error 
            % u = u(1:n-k) - phi_kk(k) * v(2:n-k+1);
            u(k+1:sampleSize) = u(k:sampleSize-1) - phi_kk(k) * v(k+1:sampleSize);
            % backward one-step prediction error 
            % v = v(2:n-k+1) - phi_kk(k) * u1(1:n-k);
            v(k+1:sampleSize) = v(k+1:sampleSize) - phi_kk(k) * u1(k:sampleSize-1);
            sigma2 = 0.5 * (1 - phi_kk(k)^2) * d/((sampleSize-k));
            d = (((1 - phi_kk(k)^2)) * d) - v(k+1)^2 - u(sampleSize)^2;
    end

    % Run Durbin-Levinson algorithm with given PACF phi_kk
    covVals = autoCorrelation(x);
    PHI(1, 1) = phi_kk(1);
    %nu(0) contains value of covariance at lag 0 i.e. covVals(1) in Matlab
    nu(1) = covVals(1) * (1 - PHI(1, 1)^2);
    for k = 2:min(order, length(covVals) - 1)
        % Take values of PACF from Burg estimates
        PHI(k, k) = phi_kk(k);
        % Calculate phi(n, 1), phi(n, 2),...,phi(n, n-1)
        PHI(k, 1:k-1) =  PHI(k - 1, 1:k-1) - PHI(k, k) * PHI(k - 1, (k-1):-1:1);
    end
    pacf = diag(PHI);
    significance  = (abs(pacf) > (1.96 / sqrt(sampleSize)));   
    phi = PHI(k, :);
    order = max(find(significance~=0));
    if (order == 0)
            phi = 1;
            pacf = 1;
    else
        phi = PHI(order,:);
        pacf = [1; pacf];
    end   
    error = sigma2;
end
