function [lambda,sigma_k] = get_lambda(sigma,N,varargin)
%GET_SIGMA Get sigma measures the sigma needed for Parzen estimate
%depending on themethod chosen
%   sigma --> sigma of random variable x
%   N --> length of the random values vector
%   multi (varargin{2}) --> dimensionality, default value is 1 (for multi = 0)
%   method (varargin{1}) --> method to calculate the lambda and sigma

if length(varargin) <= 1 
    dim = 1;
else
    dim = varargin{2};
end
if (isempty(varargin)) | (varargin{1} == 'silverman')
    sigma_k = sigma*(4*(N^-1)*((2*dim+1)^-1))^(1/(dim+4));
    lambda = 1/(2*(sigma_k^2));
end

end

