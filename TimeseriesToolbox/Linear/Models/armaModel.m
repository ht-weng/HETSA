function y = armaModel(x, thetas, phis, mu)
    if nargin() == 3
        mu = 0;
    end
    y = armaFilter([1 thetas], [1 -phis], x) + mu;
end