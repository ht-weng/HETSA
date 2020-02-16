function [pacf, significance] = partialAutoCorrelation(X)
    if(size(x, 1) > size(x, 2))
         disp('Make sure that timeseries are stored as rows');
    end
    % remove the mean value
    X = X - mean(X);
    lengthX = length(X);
    [~, ~, pacf, significance ] = estimateARbyDurbinLevinson(X, lengthX);
end