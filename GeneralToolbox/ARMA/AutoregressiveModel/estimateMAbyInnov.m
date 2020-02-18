function [theta, error] = estimateMAbyInnov(x, order)
% Estimates the parameters of the MA model using the innovation algorithm
% X - timeseries
% N - number of steps
% error - mean squared errors
% theta - parameters of the MA model
    gamma = autoCorrelation(x);
    kappa = toeplitz(gamma(1:order + 1));
    
    v(1)        = kappa(1,1); 
    THETA(1,1)  = kappa(2,1)/v(1);
    v(2)        = kappa(2,2) - THETA(1,1)^2 * v(1);
    
    for n = 2:order
            THETA(n,n) = kappa(n+1,1)/v(1);
            
            for k = 1:n-1
                    theta_kj    = THETA(k,k:-1:1);
                    theta_nj    = THETA(n,n:-1:n - k + 1);
                    THETA(n,n-k)= (kappa(n+1,k + 1) - sum(theta_kj .* theta_nj .* v(1:k)))/v(k + 1);
            end
            
            v(n+1) = kappa(n+1,n+1) - THETA(n,n:-1:1).^2 * v';
    end
    theta = THETA(end, :); 
    error = v(end);
end
