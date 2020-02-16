function x = movingaverageModel(z, mu, thetas)
    zLen = length(z);
    thetasLen = length(thetas);
    
	% pre-locate memory
    x = mu*ones(zLen, 1);
    
    for k = 1:zLen
        x(k) = x(k) + z(k);
        for l = 1:min(thetasLen, k - 1)
            x(k) = x(k) + thetas(l) * z(k - l);
        end
    end
end