function [r] = generate_gaussian_samples(N,mu,rho)
%GENERATE_INDEPENDENT This function creates N random gaussian samples
%   N indicates the number of random samples to generate
%   mu indicates the mean, if this variable is multidimensional with dimension n then n variables will be created 
%   rho is a parameter to control the correlation when mutilidemsional, if 0 then independence 
    cov = [sqrt(1-rho^2) rho; rho sqrt(1-rho^2)];
    r = mvnrnd(mu, cov, N);

end

