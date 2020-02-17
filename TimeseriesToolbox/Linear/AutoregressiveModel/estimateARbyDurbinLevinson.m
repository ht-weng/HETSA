function [phi, error, pacf, significance ]= estimateARbyDurbinLevinson(x, order)
   %  INPUT PARAMETERS:
   %  x: input time-series 
   %  order: the order of AR model to estimated
   %  significance  of estimated PACF, i.e. PHI(j, j)
   %  OUTPUT PARAMETERS:
   %  phi: the estimated coefficietns of the AR(order) model
   %  error: the mean squared error of the AR(order) model i.e. represents
   %  variance of what is left in the time-series after removing the model
   %  pacf: Partial AutoCorrelation Function, i.e.  PHI(j, j)
   %  significance : estimates significance  of each PACF value
   %  Test: page 143
   %  [phi, error, pacf, signicinance] = dl([0.17992, 0.0759, 0.04885], 2, 77)
   %  [phi, error, pacf, signicinance] = dl([1382.2, 1114.4, 591.73, 96.216], 3, 100)
   
    if(size(x, 1) > size(x, 2))
         disp('Make sure that timeseries are stored as rows');
    end
   
   sampleSize = length(x);
   gamma = autoCorrelation(x);
   order = min(order, size(gamma,2) - 1);
   PHI = zeros(order, order);
   sigma2 = zeros(order, 1);
   
   PHI(1, 1) = gamma(2)/gamma(1);
   %nu(0) contains value of covariance at lag 0 i.e. covVals(1)
   sigma2(1) = gamma(1) * (1 - PHI(1, 1)^2);
   for n = 2:min(order, length(gamma) - 1)
       % Calculate values of PACF 
       PHI(n, n) = (gamma(n + 1) - PHI(n - 1, 1:n-1) * gamma(n:-1:2)')/sigma2(n - 1);
       % Calculate phi(n, 1), phi(n, 2),...,phi(n, n-1)
       PHI(n, 1:n-1) =  PHI(n - 1, 1:n-1) - PHI(n, n) * PHI(n - 1, (n-1):-1:1);
       % Calculate prediction error of the model of order n
       sigma2(n) = sigma2(n - 1) * (1 - PHI(n, n)^2);
   end
   
   pacf = diag(PHI);
   significance  = (abs(pacf) > (1.96 / sqrt(sampleSize)));
   phi = PHI(order,:);
   order = max(find(significance~=0));
   if (order == 0)
        phi = 1;
        pacf = 1;
   else
       phi = PHI(order,:);
       pacf = [1; pacf];
   end
   error = sigma2(order);  
end

