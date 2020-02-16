function x = autoregressiveModel(z, mu, phis)
    zLen = length(z);
    phisLen = length(phis);
    
    
    % pre-locate memory
    x = mu*ones(zLen, 1);

    for k = 1:zLen
        x(k) = x(k) + z(k);
        for l = 1:min(phisLen, k - 1)
            x(k) = x(k) + phis(l) * (x(k - l) - mu);
        end
    end
end