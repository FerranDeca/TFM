function [I_2] = I2_CS(r,sigma,method,multi)
%I2_CS get an estimation of the 2-Rényi mutual information for the Chauchy-
% Schwarz method
%   r is the matrix of random values, each column a different random value
%   (Nxn)
%   sigma is the standard deviation for each column of r (nx1)
%   method must be an struct, where struct.name indicates the method to use
%   and struct.error indicates the margianl error for Cholesky
%   multi is >0 if multidimensionality is needed

    if strcmp(method.name,'kernel_matrix') 
        [mkernel,mkernel_sum,sigma_k] = prepare_estimation(r,sigma,multi);

        N = length(mkernel(:,:,1));
        n = length(r(1,:));

        L = tril(ones(N,N),-1);
        for o = 1:n
            mkernel(:,:,o) = mkernel(:,:,o).*L;
        end

        marginal = zeros(1,n);
        for o = 1:n
            marginal(o) = sum(sum((mkernel(:,:,o).^2)));
        end

        total = trace(mkernel(:,:,1)'*mkernel(:,:,2));

        I_2 = (total^2)/(prod(marginal));    
       
    elseif strcmp(method.name,'cholesky') 

        [G1,G2] = prepare_cholesky(r,sigma,method.error,multi);
        N = size(G1,2);
        L = tril(ones(N,N),-1);
        den_1= (ones(1,N)*((G1'*G1).*(G1'*G1).*L)*ones(N,1));
        den_2= (ones(1,N)*((G2'*G2).*(G2'*G2).*L)*ones(N,1));
        num = (sum(sum(((G1'*G1).*(G2'*G2))-eye(N,N))))/2;
        I_2 = ((num^2)/(den_1*den_2));
        
        
    else
        disp('Incorrect method name or empty method name input!');

    end
            
end

