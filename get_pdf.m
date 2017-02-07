function [pdf] = get_pdf(x,r,sigma)
%GET_PDF This function returns the pdf as a vector
%   x --> region where the pdf is avualuated
%   r --> random samples
%   sigma --> sigma values for the Parzen estimate

pdf = zeros(1,length(x));
N = length(r);

for i = 1:N
    kernel = generate_kernel(x,r(i),sigma);
    pdf = pdf + (1/N)*(1/(sqrt(2*pi)*sigma))*kernel;

end

