function [I_2] = I2_Parzen(r,sigma,method,multi)
%I2_Parzen get an estimation of the 2-Rényi mutual information for the method 
%   derived from the Parzen window estimate method
%   r -->   is the matrix of random values, each column a different random value
%   (Nxn)
%   sigma -->  is the standard deviation for each column of r (nx1)
%   method --> must be an struct, where struct.name indicates the method to use
%   and struct.error indicates the margianl error for Cholesky
%   multi --> is >0 if multidimensionality is needed

        [mkernel,mkernel_sum,sigma_k] = prepare_estimation(r,sigma,multi);


        for o = 1:length(r(1,:))
            mkernel(:,:,o)= mkernel(:,:,o).*(1/(sqrt(2*pi)*sigma_k(o)));
        end



        total = 0;
        for j = 1:(length(r(:,1)))-multi+1
            num = sum(prod(mkernel(j,:,:),3))/(length(r(:,1))-multi+1);
            den = prod(sum(mkernel(j,:,:)),3)*((length(r(:,1))-multi+1)^(-length(r(1,:))));
            total = num/den + total;
        end
        I_2 = log10(total/((length(r(:,1)))-multi+1)); 


end

