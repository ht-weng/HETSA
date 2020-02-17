function y = armaFilter(thetas, phis, x)
    % direct form II transposed
    if(size(x, 1) > size(x, 2))
        disp('Make sure that timeseries are stored as rows');
    end
    
    phis_len = length(phis);
    thetas_len = length(thetas);
    inLen = size(x,2);
	y = zeros(size(x));
%     initialValue = x(1);
%     x = x - initialValue;
    for k = 1:size(x,1)
        tmp = 0;
        for i = 1:inLen
            y(k, i) = 0;
            % Moving average
            for j = 1:thetas_len
                if(i - j < 0) 
                    continue;
                end
                tmp = tmp + thetas(j) * x(k, i - j + 1);
            end
            % Autoregressive part
            for j = 2:phis_len
                if(i - j < 0) 
                    continue;
                end
                tmp = tmp - phis(j) * y(k, i - j + 1);
            end

            tmp = tmp / phis(1);
            y(k, i) = tmp;
            tmp = 0;
        end
    end
%    y = y + initialValue;
end