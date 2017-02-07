function [mkernel,mkernel_sum,sigma_k] = prepare_estimation(r,sigma,multi)
%PREPARE_ESTIMATION This function measures the estimation of the pdf, the
%estimation of the Renyi entropy and the Kernel matrix
%   r --> matrix fo the two random variables
%   sigma --> sigma of the data
%   multi --> dimensionality indicator, for point to point this value is 1

    n = length(r(1,:));
    x = -10:0.05:10;
    N = length(r(:,1));
    
    if multi < 2

        pdf = zeros(n,length(x));
        for o = 1:n
            pdf(o,:) = get_pdf(x,r(:,o),sigma(o));
        end

        mkernel_sum = zeros(n,1);
        mkernel = zeros(N,N,n);
        lambda = zeros(n,1);
        sigma_k = zeros(n,1);
        for o = 1:n
            [lambda(o),sigma_k(o)] = get_lambda(sigma(o),N,'silverman',multi);
        end
        for o = 1:n
            [mkernel(:,:,o),mkernel_sum(o)] = get_mkernel(lambda(o),r(:,o),[1 N;1 N],multi);
        end

    else

        pdf = zeros(n,length(x));
        for o = 1:n
            pdf(o,:) = get_pdf(x,r(:,o),sigma(o));
        end
        mkernel_sum = zeros(n,1);
        mkernel = zeros(N-multi+1,N-multi+1,n);
        lambda = zeros(n,1);
        sigma_k = zeros(n,1);        

        for o = 1:n
            [lambda(o),sigma_k(o)] = get_lambda(sigma(o),N-multi+1,'silverman',multi);
        end

        for o = 1:n
            [mkernel(:,:,o),mkernel_sum(o)] = get_mkernel(lambda(o),r(:,o),[1 N;1 N],multi);
        end
        
    end

end

