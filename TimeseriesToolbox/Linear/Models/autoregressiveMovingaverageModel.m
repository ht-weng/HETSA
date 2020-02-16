function x = autoregressiveMovingaverageModel(z, mu, thetas, phis)
    zLen = length(z);
    thetasLen = length(thetas);
    phisLen = length(phis);
    
	% pre-locate memory
    x = mu*ones(zLen, 1);
    
    for k = 1:zLen
        % MA part
        x(k) = x(k) + z(k);
        for l = 1:min(thetasLen, k - 1)
            x(k) = x(k) + thetas(l) * z(k - l);
        end
        % AR part
        for l = 1:min(phisLen, k - 1)
            x(k) = x(k) + phis(l) * (x(k - l) - mu);
        end
    end
end