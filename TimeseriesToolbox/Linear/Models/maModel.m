function y = maModel(x, thetas, mu)
    if nargin() == 2
        mu = 0;
    end
    y = armaFilter([1 thetas], 1, x) + mu;
end