function y = arModel(x, phis, mu)
    if nargin() == 2
        mu = 0;
    end
    y = armaFilter(1, [1 -phis], x) + mu;
end