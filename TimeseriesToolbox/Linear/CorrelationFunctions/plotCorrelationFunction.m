function plotCorrelationFunction(x, N)
    if(size(x, 1) > size(x, 2))
         disp('Make sure that timeseries are stored as rows');
    end
    lengthX = length(x);
    if(nargin == 1)
        N = round(0.1*lengthX);
    end
    
    figure, stem(0:(N-1), x(1:N)); % plot CF with confidence interval
    line([0, N], [1.96/sqrt(lengthX), 1.96/sqrt(lengthX)], 'color', 'r');
    line([0, N], [-1.96/sqrt(lengthX), -1.96/sqrt(lengthX)], 'color', 'r');
end