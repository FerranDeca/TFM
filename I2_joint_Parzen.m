function [I_2] = I2_joint_Parzen(r,sigma,method,multi)
%I2_Parzen get an estimation of the 2-Rényi mutual information for the method 
%   derived from lack of sub-additivity 
%   r -->   is the matrix of random values, each column a different random value
%   (Nxn)
%   sigma -->  is the standard deviation for each column of r (nx1)
%   method --> must be an struct, where struct.name indicates the method to use
%   and struct.error indicates the margianl error for Cholesky
%   multi --> is >0 if multidimensionality is needed

    N = length(r);
    n= size(r,2);

    if strcmp(method.name,'kernel_matrix')
        [mkernel,mkernel_sum,sigma_k] = prepare_estimation(r,sigma,multi);

%         for o = 1:n
%             for i = 1:N
%                 mkernel(i,i,o) = 0;
%             end
%         end
%         
%         M = N-multi+1;
%         num = sum(sum(prod(mkernel,3)));
%         total = num/(M*(M-1));
%         total2 = ((1/(M*(M-1)))^2)*prod(mkernel_sum);
%         I_2 = abs(total-total2);
%         stem((total/(prod(mkernel_sum))));
%         hold on;
        
        
        num = sum(sum(prod(mkernel,3)));
        total = num*((N-multi+1)^2);
        I_2 = abs((total/(prod(mkernel_sum)))-1);
%         stem((total/(prod(mkernel_sum))));
%         hold on;

    elseif strcmp(method.name,'cholesky')
        

        [G1,G2] = prepare_cholesky(r,sigma,method.error,multi);
        total= norm(G2*G1','fro')^2;
        num_1= norm(ones(1,N-multi+1)*G1','fro')^2;
        num_2= norm(ones(1,N-multi+1)*G2','fro')^2;
        I_2 = abs((((N-multi+1)^2)*total/(num_1*num_2))-1);


    else
        disp('Incorrect method name or empty method name input!');
    end
    
    
end

