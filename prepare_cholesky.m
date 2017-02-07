function [G1,G2] = prepare_cholesky(r,sigma,error,multi)
%PREPARE_CHOLESKY This function returns the cholesky decomposition for a
%given matrix of random values r (Nx2)
%   r --> matrix fo the two random variables
%   sigma --> sigma of the data
%   multi --> dimensionality indicator, for point to point this value is 0

    n = length(r(1,:));
    N = length(r(:,1));
    lambda = zeros(n,1);
    sigma_k = zeros(n,1);
    for o = 1:n
        [lambda(o),sigma_k(o)] = get_lambda(sigma(o),N,'silverman',multi);
    end
    
    if multi < 2
        for o = 1:n
            G = zeros(N,N);
            [G(:,1),m_sum]=get_mkernel(lambda(o),r(:,o),[1 N; 1 1],multi);
            j = 1;
            D = zeros(N,1);
            P = (1:1:N)';
            for i = 1:N
                if i == 1
                    D(i:N) = ones(N,1);
                else
                    D(i:N) = ones(N-i+1,1)-(G(i:N,1:i-1).*G(i:N,1:i-1))*ones(i-1,1);
                end
                suma = sum(D(i:end));
                if suma < error
                    break;
                end
                [val,idx]= max(D(i:end));
                if i > 1
                    j=idx+i-1;
                end


                temp = P(i); 
                P(i) = P(j);
                P(j) = temp;

                temp_2 = G(i,1:i-1); 
                G(i,1:i-1) = G(j,1:i-1);
                G(j,1:i-1) = temp_2;

                G(i,i)= sqrt(D(j));

                [mkernel_pivot,m_sum]=get_mkernel(lambda(o),r(:,o),[i+1 N; i i],multi,P);
                G(i+1:N,i)=(mkernel_pivot-G(i+1:N,1:i-1)*(G(i,1:i-1))')/(G(i,i));
            end
            [sortedA,IX] = sort(P);
            sortedB = G(IX,:);
            G_1 = sortedB(:,1:i-1);
            if o == 1
                G1=G_1';
            else
                G2=G_1';
            end
        end
        
    else
        
        for o = 1:n
            N = length(r(:,1));
            N2 = N-multi+1;
            G = zeros(N2,N2);
            P = (1:1:N2)';
            [G(:,1),m_sum]=get_mkernel(lambda(o),r(:,o),[1 N; 1 1],multi);
            j = 1;
            D = zeros(N2,1);

            for i = 1:N2
                if i == 1
                    D(i:N2) = ones(N2,1);
                else
                    D(i:N2) = ones(N2-i+1,1)-(G(i:N2,1:i-1).*G(i:N2,1:i-1))*ones(i-1,1);
                end
                suma = sum(D(i:end));
                if suma < error
                    break;
                end
                [val,idx]= max(D(i:end));
                if i > 1
                    j=idx+i-1;
                end

                P([i j])=P([j i]);

                temp_2 = G(i,1:i-1); 
                G(i,1:i-1) = G(j,1:i-1);
                G(j,1:i-1) = temp_2;

                G(i,i)= sqrt(D(j));

%                 [mkernel_pivot,m_sum]=get_mkernel(lambda(o),r(:,o),[i+1 N; 1 1],multi,P);
                [mkernel_pivot,m_sum]=get_mkernel(lambda(o),r(:,o),[i+1 N; i i],multi,P);
                
                G(i+1:N2,i)=(mkernel_pivot-G(i+1:N2,1:i-1)*(G(i,1:i-1))')/(G(i,i));
            end
            [sortedA,IX] = sort(P);
            sortedB = G(IX,:);
            G2 = sortedB(:,1:i-1);

            if o == 1
                G1=G2';
            else
                G2=G2';
            end
        end
        
        
        
        
    end
        
        




end

