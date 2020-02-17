function y = armaFilter2(thetas, phis, x)
    % direct form II transposed
    xLen = size(x, 2);
    thetasLen = length(thetas);
    phisLen = length(phis);
    
	% pre-locate memory
    y = zeros(size(x));
    
    for k = 1:xLen
        for m = 1:size(x,1)
            % MA part
            y(m, k) = sum(thetas(1:min(k, thetasLen)) .* x(m, k:-1:max(k-thetasLen + 1, 1)));
            % AR part
            y(m, k) = -y(m, k);
            y(m, k) = -sum(phis(1:min(k, phisLen)) .* y(m, k:-1:max(k-phisLen + 1, 1)));
        end
    end
end