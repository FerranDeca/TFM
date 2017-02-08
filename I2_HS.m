function [I_2] = I2_HS(r,sigma,method,multi)

    N = length(r(:,1));
    n = length(r(1,:));
    
    if strcmp(method.name,'kernel_matrix') 

        [mkernel,mkernel_sum,sigma_k] = prepare_estimation(r,sigma,multi);

    %     
    %     L = triu(ones(N,N),-1);
    %     mkernel2 = mkernel;
    %     for o = 1:n
    %         mkernel2(:,:,o) = mkernel(:,:,o).*L;
    %     end

        for o = 1:n
            for i = 1:N-multi+1
                mkernel(i,i,o) = 0;
            end
        end


%         P = eye(N,N) - (1/N)*ones(N,1)*ones(N,1)';
%         numerador= trace((P*mkernel(:,:,1)*P*mkernel(:,:,2)*P)'*(P*mkernel(:,:,1)*P*mkernel(:,:,2)*P));
%         I_2 = (1/N)*sqrt(numerador);


        numerador= trace(mkernel(:,:,1)*mkernel(:,:,2));
        I_2 = (1/(N-multi+1))*sqrt(numerador);

    elseif strcmp(method.name,'cholesky') 

        [G1,G2] = prepare_cholesky(r,sigma,method.error,multi);
        numerador = trace((G1'*G1)*(G2'*G2))-N-multi+1;
        I_2 = (1/(N-multi+1))*sqrt(numerador);
        

    else
        disp('Incorrect method name or empty method name input!');

    end


end

