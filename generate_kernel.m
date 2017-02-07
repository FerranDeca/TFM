function kernel = generate_kernel(x,x_,sigma)
%GENERATE_KERNEL Function that returns a Gaussian Kernel with standard
%deviation based on the sigma parameter and centered at x_
%   x -->               region where the kernel will be evaluated
%   x_ -->              bias of the kernel position
%   sigma -->           standard deviation
%   output kernel -->   vector containing the kernel generated with the
%                       parameters introduced
%   If x is not a vector but a point, then the function returns the
%   distance between two point in a Gaussian kernel

    lambda = 1/(2*(sigma^2));
    kernel = exp(-lambda*(x-x_).^2);

end

